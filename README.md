# OpenFOAM library for physiological flow

A library of functions and tools to solve for physiological flow in OpenFOAM. The following is included

* Editted pimpleFoam solver that can compute 3 element Windkessel model boundary conditions. 

## Getting Started

### Prerequisites

These files have been developed with OpenFOAM 4.0. 
This can be downloaded from http://www.openfoam.com/

### Installing

Once you have the software installed, download the folder that you want to use.

If it is a solver, copy it into an apps folder. 
If you want, change the name of the folder and its name in Make/files script. 

Then inside the folder, with OpenFOAM loaded, perform the following commands:

wclean
wmake

You can check if is has compiled correctly by using

ls $FOAM_USER_APPBIN

If the solver is listed, then it is ready to use. 

### Deployment

The sample case file shows how it can be implemented. 




