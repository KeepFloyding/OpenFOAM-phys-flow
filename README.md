# OpenFOAM library for physiological flow

A library of solvers and boundary conditions to solver for physiological flow in OpenFOAM.

The following solvers are included:

* windkesselSolver (editted pimpleFoam solver that can implement 3 element Windkessel boundary conditions)

The following boundary conditions are included:

* parabolicVelocity (inlet condition that specifies a parabolic velocity distribution)

## Getting Started

### Prerequisites

These files have been developed with OpenFOAM 4.0. 
This can be downloaded from http://www.openfoam.com/

Alternatively you may be able to load the module with 

module load openfoam/4.0

Set up an openfoam directory if you don't have one with

mkdir -p $FOAM_RUN

cd $FOAM_RUN

cd ..

FOAM_DEV=$PWD

mkdir -p $FOAM_DEV

### Installing

For solvers:

cd $FOAM_DEV

git clone https://github.com/KeepFloyding/OpenFOAM/tree/master/[SOLVER_NAME]

cd [SOLVER_NAME]/[SOLVER_VERSION]

wclean

wmake

// Check if the solver is available

ls $FOAM_USER_APPBIN

For boundary conditions:

cd $FOAM_DEV

git clone https://github.com/KeepFloyding/OpenFOAM/tree/master/boundaryConditions

cd boundaryConditions

wclean

wmake libso

// Check if the boundary conditions are available

ls $FOAM_USER_LIBBIN

### Deployment

If compilation is successful, then the solver and boundary conditions can be used.
The sample case file shows how they can be implemented. 

For solvers:
Create case file similar to sampleCaseFile and launch by typing [SOLVER_NAME]

For boundary conditions:
Edit the system/controlDict file and add at the end

libs ("[LIB_NAME]");

Implementation can be seen in the case file








