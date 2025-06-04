/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

  License
  This file is part of ITHACA-FV

  ITHACA-FV is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ITHACA-FV is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Publ  ic License
  along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "unsteadyNS.H"

/// \file
/// Source file of the unsteadyNS class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
//for initializing
unsteadyNS::unsteadyNS() {}

// Construct from zero
unsteadyNS::unsteadyNS(int argc, char* argv[])
    :
    UnsteadyProblem()  //setTimes() and checkWrite()
{
    _args = autoPtr<argList>  //create a new argList object, have _args() to query the parsed runtime arugments
            (
                new argList(argc, argv)  //construct an argList instance
            );

    if (!_args->checkRootCase())  //if checkRootCase() returns false,then the directory setup isn't a proper case
    {
        Foam::FatalError.exit();  //then FatalError.exit() is called
    }

    argList& args = _args();  //now we can use args. rahter than _args()->
#include "createTime.H"
#include "createMesh.H"
    _pimple = autoPtr<pimpleControl>  //call the pimpleControl constructor, passing in mesh
              (
                  new pimpleControl
                  (
                      mesh
                  )
              );
    ITHACAdict = new IOdictionary  //construct an IOdictionary from IOobject
    (
        IOobject
        (
            "ITHACAdict",  //the filename to look for
            runTime.system(),  //point to system folder
            mesh,  //associate it with mesh database so that any mesh-dependent lookups can be resolved correctly
            IOobject::MUST_READ,  //use MUST_READ to ensure the solver will abort if ITHACAdict not present
            IOobject::NO_WRITE  //dictionary is read-only
        )
    );
#include "createFields.H"
#include "createFvOptions.H"
    para = ITHACAparameters::getInstance(mesh, runTime);  //ITHACAparameters is the class designed to hold all user-defined parameters
    //if ITHACAparameters object doesn't exist in memory, then getInstance(mesh, runTime) will construct one by reading a dictionary
    //using mesh and runTime context
    method = ITHACAdict->lookupOrDefault<word>("method", "supremizer");
    //find a keyword called method in ITHACAdict, if exist, return a string, otherwise return "supremizer"
    M_Assert(method == "supremizer"
             || method == "PPE",
                       "The method must be set to supremizer or PPE in ITHACAdict");
    //if method is neither supremizer nor PPE, then aborts and prints "The method must be set..."; otherwise, execution continues normally
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty",
             "The BC method must be set to lift or penalty in ITHACAdict");
    timedepbcMethod = ITHACAdict->lookupOrDefault<word>("timedepbcMethod", "no");
    M_Assert(timedepbcMethod == "yes" || timedepbcMethod == "no",
             "The BC method can be set to yes or no");
    timeDerivativeSchemeOrder =
        ITHACAdict->lookupOrDefault<word>("timeDerivativeSchemeOrder", "second");
    M_Assert(timeDerivativeSchemeOrder == "first"
             || timeDerivativeSchemeOrder == "second",
                                          "The time derivative approximation must be set to either first or second order scheme in ITHACAdict");
    offline = ITHACAutilities::check_off();  //Check if the offline data folder "./ITHACAoutput/Offline" exists
    podex = ITHACAutilities::check_pod();  //Check if the POD data folder "./ITHACAoutput/POD" exists
    supex = ITHACAutilities::check_sup();  //Check if the supremizer folder exists
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void unsteadyNS::truthSolve(List<scalar> mu_now, fileName folder)  //mu_now is a vector of parameter values like Reynolds number, boundary condition coefficients, etc.
{
    Time& runTime = _runTime();  //manages current time, delta t and writing directories
    surfaceScalarField& phi = _phi();  //phi represents the volumetric flux on each face of control volumes
    fvMesh& mesh = _mesh();  //use mesh to represent _mesh()
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField& p = _p();
    volVectorField& U = _U();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();  //instant is a scalar for physical time value, and runTime.times() reads and stores these instants
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;  //nextWrite tracks the absolute time when we should export next field

    // Set time-dependent velocity BCs for initial condition
    if (timedepbcMethod == "yes")  //set timedepbcMethod in ITHACAdict
    {
        for (label i = 0; i < inletPatch.rows(); i++)  //number ofinlet pathces we defined
        {
            Vector<double> inl(0, 0, 0);

            for (label j = 0; j < inl.size(); j++)
            {
                inl[j] = timeBCoff(i * inl.size() + j, 0);  //compute which row holds the j-th component of patch i, pick column 0
            }

            assignBC(U, inletPatch(i, 0), inl);  //wirte inl into BC for field U on that patch(i,0)
        }
    }

    // Export and store the initial conditions for velocity and pressure
    ITHACAstream::exportSolution(U, name(counter), folder);  //write the field U or P into a new time-snapshot directory
    ITHACAstream::exportSolution(p, name(counter), folder);  //counter=0,1,2,... to ensure U and P at the currrent physical time are saved in a uniquely named floder
    std::ofstream of(folder + name(counter) + "/" +
                     runTime.timeName());  //example: ./ITHACAoutput/Offline/0/0/0
    Ufield.append(U.clone());  //U.clone() makes a deep copy of volVecotrField U at this moment
    Pfield.append(p.clone());
    counter++;
    nextWrite += writeEvery;  //writeEvery is the physical time of snapshots

    // Start the time loop
    while (runTime.run())  //read controls, update delta t, solve momentum+pressure, write snapshots until endTime
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;  //timeIndex() plus 1, timeValue() += deltaT()
        Info << "Time = " << runTime.timeName() << nl << endl;  //print physical time step

        // Set time-dependent velocity BCs
        if (timedepbcMethod == "yes")
        {
            for (label i = 0; i < inletPatch.rows(); i++)
            {
                Vector<double> inl(0, 0, 0);

                for (label j = 0; j < inl.size(); j++)
                {
                    inl[j] = timeBCoff(i * inl.size() + j, counter2);  //counter2 is the column to pick up
                }

                assignBC(U, inletPatch(i, 0), inl);
            }

            counter2 ++;  //for next iteration we pull the next-in-time column of timeBCoff
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())  //pimple.loop()=TRUE if the outer PIMPLE iteration count hasn't reached maximum
        {
#include "UEqn.H"  //solve for provisional velocity

            // --- Pressure corrector loop, to enforce continuity
            while (pimple.correct())  //pimple.correct()=TRUE if inner loop converges
            {
#include "pEqn.H"  //build and solve a pressure Poisson-type eqation
            }

            if (pimple.turbCorr())  //update turbulence
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"  //processor time during last time-step solve
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"  //real time since last time-step began
             << nl << endl;

        if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            Ufield.append(U.clone());
            Pfield.append(p.clone());
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            //mu_samples each row records one snapshot's "physical time"(first column) and all parameter values(next columns)
            //conservativeResize adds one more row to bottom
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());
            //the first entry in new row(column 0) is current time, convert from string to double
            for (label i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];  //have one row forevery time the solver wrote a snapshot
            }
        }
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)  //never set up any parameter
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == counter * mu.cols())  ..successfully collect exactly one row for every snapshot
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   folder);
    }
}
