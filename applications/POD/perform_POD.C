/*---------------------------------------------------------------------------*\
Copyright (C) 2017 by the ITHACA-FV authors

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

    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

Class
    POD

Description
    Application to perform a POD decomposition of a general field on a given case

SourceFiles
    perform_POD.C

\*---------------------------------------------------------------------------*/

/// \file
/// \brief Application to perform POD on an already run case
/// \details In order to use this file one needs to prepare a ITHACAPODdict file, in order to 
/// check the syntax one needs to check the \ref ITHACAPODdict file.

/// \file ITHACAPODdict
/// \brief Example of a ITHACAPODdict file


#include "fvCFD.H"  //core finite-volume classes
#include "IOmanip.H"  //routines that build the snapshot correlation matrix and solve for eigenmodes
#include "IFstream.H"  
#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "volFieldsFwd.H"
#include "forces.H"
#include "forceCoeffs.H"
#include "volFields.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include "ITHACAPOD.H"
#include "ITHACAparameters.H"  //singeleton holding things like boundary-condition correction flags
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])  //argc is the number of arguments, argv is the array of arguments
{

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

    PtrList<volVectorField> Vfield;  //a list of pinters to objects of type <volVectorField>
    PtrList<volScalarField> Sfield;  //Sfield will hold in-memory clones of each snapshot we read, so does Vfield
    PtrList<volVectorField> Vmodes;  //Vmodes will store the resulting POD modes returned by ITHACAPOD::getModes()
    PtrList<volScalarField> Smodes;  //same as Vmodes

    ITHACAparameters* para = ITHACAparameters::getInstance(mesh,
                             runTime);
    //pointer para ponits to the type of data which is ITHACAparameters.
    //for right side, get getInstance from ITHACAparameters, passing in mesh and the simulation runTime to read relevant information

    bool pod_exist;  //declares a plain C++ Boolean variable
    struct stat sb;  //stat is a POSIX system call that fills in a struct stat with information about a file or directory
                     //sb is the variable that will hold all metadata

    if (stat("./ITHACAoutput/POD", &sb) == 0 && S_ISDIR(sb.st_mode))
    //modifying sb inside function stat will at the same time modify sb in the function that called stat
    {
        pod_exist = true;
    }
    else
    {
        pod_exist = false;
        Info << "POD don't exist, performing a POD decomposition" << endl;
        mkDir("./ITHACAoutput/POD");  //make sure there's a POD folder under ITHACAoutput/ before starting dumping modes
        system("ln -s ../../constant ./ITHACAoutput/POD/constant");  //execute arbitrary shell commands directly from the code
        system("ln -s ../../0 ./ITHACAoutput/POD/0");
        system("ln -s ../../system ./ITHACAoutput/POD/system");
    }
    if(pod_exist == 1)
    {
        Info << "The POD has already been performed, please delete the ITHACAoutput folder and try again." << endl;
        exit(0);
    }

    // Initialize Variables
    label nSnapshots;  //label is a typedef for an integer type int, nSnapshots is total numver of snapshots(time-step)
    label startTime;  // index of the first snapshot
    label endTime;

    //OpenFOAM Read FORCESdict
    IOdictionary ITHACAPODdict
    (
        IOobject
        (
            "ITHACAPODdict",  //sepcifies the file name in the case's system/ directory
            runTime.system(),  //return the path to system/ in the case, that is looking for ./system/ITHACAPODdict
            mesh,
            IOobject::MUST_READ,  //if not find ITHACAPODdict, OpenFOAM will immediately throw an error and abort
            IOobject::NO_WRITE  //no modifying
        )
    );

    // Get times list from the case folder
    instantList Times = runTime.times();
    //instantList is a typedef for List<scalar/double>

    // Read Initial and last time from the POD dictionary
    const entry* existnsnap = ITHACAPODdict.lookupEntryPtr("Nsnapshots", false, true); //number of snapshots
    const entry* existLT = ITHACAPODdict.lookupEntryPtr("FinalTime", false, true);  //set the time window
    //lookupEntryPtr returns a raw pointer to the internal ectry object for the key(Nsnapshots & FinalTime), or nullptr if not found
    //first Bool is flase if nothing missing
    //second Bool is for deep searching

    // Initiate variable from PODSolverDict
    if ((existnsnap) && (existLT))
    {
        Info << "Error you cannot define LatestTime and NSnapShots together" << endl;
        abort();
    }
    else if (existnsnap)  //only when dictionary specified a fixed number of snapshots
    {
        scalar InitialTime = ITHACAPODdict.lookupOrDefault<scalar>("InitialTime", 0);  //give start time. If no InitilaTime, it goes to 0
        nSnapshots = readScalar(ITHACAPODdict.lookup("Nsnapshots"));
        startTime = Time::findClosestTimeIndex(runTime.times(), InitialTime);  //find the closest time to InitialTime
        nSnapshots = min(nSnapshots , Times.size() - startTime);  //make sure not run past the end of time list
        endTime = startTime + nSnapshots - 1;
        Info << nSnapshots << endl;
    }
    else  //set up snapshot window only by InitialTime and FinalTime
    {
        scalar InitialTime = ITHACAPODdict.lookupOrDefault<scalar>("InitialTime", 0);
        scalar FinalTime = ITHACAPODdict.lookupOrDefault<scalar>("FinalTime", 100000000000000);
        endTime = Time::findClosestTimeIndex(runTime.times(), FinalTime); 
        startTime = Time::findClosestTimeIndex(runTime.times(), InitialTime);
        nSnapshots = endTime - startTime + 1;
        if (InitialTime > FinalTime)
        {
            Info << "FinalTime cannot be smaller than the InitialTime check your ITHACAPODdict file\n" << endl;
            abort();
        }
    }
    // Print out some Infos
    Info << "startTime: " << startTime << "\n" << "endTime: " << endTime << "\n" << "nSnapshots: " << nSnapshots << "\n" << endl;

    // Set the initial time
    runTime.setTime(Times[startTime], startTime);

    wordList fieldlist
    (
        ITHACAPODdict.lookup("fields")
    );

    //word Name = ITHACAPODdict.lookup("fieldName");
    //word type = ITHACAPODdict.lookup("type");

    if (startTime == endTime)
    {
        Info << "The case has no snapshots to process, exiting the code" << endl;
        exit(0);
    }

    for (label k = 0; k < fieldlist.size(); k++)  //k runs once per entry in fields list
    {
        dictionary& subDict = ITHACAPODdict.subDict(fieldlist[k]);  //sub-dictionary for the field name, like "U_pod"
        scalar nmodes = readScalar(subDict.lookup("nmodes"));
        word field_name(subDict.lookup("field_name"));
        word field_type(subDict.lookup("field_type"));

        for (label i = startTime; i < endTime + 1; i++)  //loop over every snapshot index in chosen time window
        {
            Info << "Reading snapshot " << i << " for field " << field_name << endl;
            //let runTime an mesh use the i-th time directory
            runTime.setTime(Times[i], i);
            mesh.readUpdate();

            if (field_type == "vector")
            {
                //read U (vector field) from disk at that time
                volVectorField vector_field
                (
                    IOobject
                    (
                        field_name,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                );
                //clone it into snapshot list
                Vfield.append(vector_field.clone());
            }

            if (field_type == "scalar")
            {
                volScalarField scalar_field
                (
                    IOobject
                    (
                        field_name,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                );
                Sfield.append(scalar_field.clone());
            }
        }

        if (field_type == "vector")
        {
            ITHACAPOD::getModes(Vfield, Vmodes, field_name, 0, 0, 0, nmodes, para->correctBC);
            //keep at zero to select the default behavior in ITHACApod
            //para->correctBC ensures the POD basis respects no-slip or other essential boundary conditions correctly
        }
        if (field_type == "scalar")
        {
            ITHACAPOD::getModes(Sfield, Smodes, field_name, 0, 0, 0, nmodes, para->correctBC);
        }
        //for getModes, it build correlation matrix using the snapshots in snapshotList
        //secondly, it finds eigenvectors and eigenvalues to solve eigneproblem
        //thirdly, it forms physical modes and sotre the first nmodes into modesList
        //lastly, dumps each mode to disk under POD/ with name fieldName

        Vfield.clear();
        Sfield.clear();
        //delete objecs and reset list size to zero, to save place for memory
    }
    Info << endl;
    Info << "End\n" << endl;
    return 0;
}


