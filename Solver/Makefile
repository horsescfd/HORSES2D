#
#///////////////////////////////////////////////////////////////////////////////////////////////////////
#
#    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
#    Copyright (C) 2017  Juan Manzanero Torrico (juan.manzanero@upm.es)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#////////////////////////////////////////////////////////////////////////////////////////////////////////
#
SYSTEM ?= LINUX
COMPILER ?=gfortran
COMM ?= SERIAL
INC_DIR=./include
BUILD_DIR=./build
BIN_DIR=./bin
SRC_DIR=./src
LIB_DIR=./lib
EXECUTABLE = HORSES2D.NS
MODE ?=DEBUG
BLAS_LIB = 
LAPACK_LIB = 
LIBS=
GLOB_MACROS := -D_DIMENSIONLESS_TAU
CURRENT_PATH := ${CURDIR}

# --------------- MACROS ---------------------------------------------------

#GLOB_MACROS +=-D_USE_LAPACK

# --------------- DEPENDENCIES ---------------------------------------------------


Physics: PhysicsNS Headers
Stopwatch: SMConstants
PhysicsNS: Setup SMConstants Headers MatrixOperations
VariableConversion: Physics
RiemannSolvers: Physics SMConstants InviscidFluxes
InviscidFluxes: Physics SMConstants
IntegerArrayLinkedList: Qsort
ViscousFluxes: Physics SMConstants
DGSEM_class: Physics SMConstants NodesAndWeightsClass Setup QuadMesh MeshFile DGSpatialDiscretizationMethods DGTimeIntegrator Storage DGBoundaryConditions Plotter Stopwatch
Plotter: SMConstants
Tecplot: Plotter
NodesAndWeightsClass: SMConstants InterpolationAndDerivatives LegendreAlgorithms MatrixOperations Physics
InterpolationAndDerivatives: SMConstants Utilities
Utilities: RealDataLinkedList SMConstants
RealDataLinkedList: SMConstants
LegendreAlgorithms: SMConstants
QuadMesh: NodeClass QuadElement MeshFile Setup InitialConditions Storage DGBoundaryConditions Utilities IntegerArrayLinkedList
MeshFile: Setup SMConstants NetCDF_Interface Headers IntegerArrayLinkedList
QuadElement: NodeClass SMConstants Physics NodesAndWeightsClass Storage InterpolationAndDerivatives MatrixOperations MetricMatrix
NodeClass: SMConstants Physics
MetricMatrix: SMConstants
Checks: DGSEM_class Physics SMConstants NodesAndWeightsClass Setup QuadMesh MeshFile DGSpatialDiscretizationMethods Storage DGBoundaryConditions
Headers: SMConstants
NetCDF_Interface: SMConstants
ProblemFile: SMConstants Setup Physics DGSEM_class
InitialConditions: SMConstants Physics
Setup: SMConstants Paramfile
DGWeakIntegrals: SMConstants Physics QuadElement MatrixOperations
DGSpatialDiscretizationMethods: SMConstants Physics QuadMesh DGViscousMethods MatrixOperations DGInviscidMethods DGWeakIntegrals DGArtificialDissipation
DGArtificialDissipation: SMConstants QuadMesh QuadElement Setup
DGViscousMethods: SMConstants Physics Setup QuadElement MatrixOperations NodesAndWeightsClass QuadMesh Headers DGWeakIntegrals
DGInviscidMethods: SMConstants Physics Setup QuadElement MatrixOperations NodesAndWeightsClass
DGViscous_BR1: DGViscousMethods
DGViscous_IP: DGViscousMethods
DGBoundaryConditions: SMConstants  Setup Physics QuadElement Paramfile
DirichletBC: DGBoundaryConditions
EulerWallBC: DGBoundaryConditions
FarfieldBC: DGBoundaryConditions
PeriodicBC: DGBoundaryConditions
PressureInletBC: DGBoundaryConditions
PressureOutletBC: DGBoundaryConditions
RiemannBC: DGBoundaryConditions
ViscousWallBC: DGBoundaryConditions
FileWriting: SMConstants Physics QuadMesh Setup
Paramfile: SMConstants
DGTimeIntegrator: SMConstants QuadMesh DGSpatialDiscretizationMethods Setup FileWriting Plotter Monitors Stopwatch
Monitors: SMConstants Setup Paramfile Physics QuadMesh
Storage: SMConstants
main: SMConstants NodesAndWeightsClass Physics MeshFile DGSEM_class Checks


# --------------- DO NOT EDIT ----------------------------------------------------

F90_REMOVE=ProblemFile.f90 Checks.f90
F90_INCLUDE=-I$(SRC_DIR)/Mesh/ -I$(SRC_DIR)/Generic/

OBJS_F90:= $(notdir $(wildcard $(SRC_DIR)/**/*.f90)) $(notdir $(wildcard $(SRC_DIR)/*.f90))
OBJS_F90:= $(filter-out $(F90_REMOVE),$(OBJS_F90))
OBJS= $(patsubst %.f90,%,$(OBJS_F90))
OBJS_EXEC = $(foreach obj,$(OBJS),$(BUILD_DIR)/$(obj).o) 


include make.inc

# -------------- CHOOSE SYSTEM -----------------------
ifeq ($(SYSTEM),MACOSX)
   DYNLIB_FLAG=-fPIC -shared -install_name @rpath/libproblemfile.dylib -o ./lib/libproblemfile.dylib
   PROBLEM_LIB=./lib/libproblemfile.dylib -Wl,-rpath,./SETUP -Wl,-rpath,$(CURRENT_PATH)/lib
else ifeq ($(SYSTEM),LINUX)
   DYNLIB_FLAG=-fPIC -shared -o ./lib/libproblemfile.so
   PROBLEM_LIB=-L./lib -lproblemfile -Wl,-rpath=./SETUP -Wl,-rpath=$(CURRENT_PATH)/lib
endif

# -------------- CHOOSE COMPILER ------------------
ifeq ($(COMPILER),gfortran)
   MOD_FLAG :=-J
   FC=gfortran
   ifeq ($(MODE),DEBUG)
      FFLAGS = -ffree-line-length-0 -cpp -O0 -g -fcheck=all -Wno-unused-variable -fbacktrace -fbounds-check -finit-real=snan -Wall -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries -fbackslash -fcheck=pointer
      CHECKFFLAGS = -ffree-line-length-0 -cpp -O0 -g -fcheck=all -Wno-unused-variable -fbacktrace -fbounds-check -finit-real=snan -Wall -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries -fbackslash -fcheck=pointer
   else 
      FFLAGS = -cpp -ffree-line-length-0 -O3 -ftree-vectorize  -ftree-vectorizer-verbose=0 -fbackslash
      CHECKFFLAGS = -cpp -ffree-line-length-0 -O3 -ftree-vectorize  -ftree-vectorizer-verbose=0 -fbackslash
   endif
   LIBS += -L$(BLAS_LIB) -L$(LAPACK_LIB) -llapack -lblas -L$(NETCDF_LIB) -lnetcdff $(PROBLEM_LIB)
else ifeq ($(COMPILER),ifort)
   MOD_FLAG := -module 
   FC=ifort
   ifeq ($(MODE),DEBUG)
     FFLAGS = -cpp -O0 -g -warn all -fpscomp logicals -mkl -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces -ftrace=full -assume ieee_fpe_flags -ftz -fp-model precise -fp-speculation=off -assume bscc
     CHECKFFLAGS = -cpp -O0 -g -warn all -fpscomp logicals -mkl -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces -ftrace=full -assume ieee_fpe_flags -ftz -fp-model precise -fp-speculation=off -assume bscc
   else ifeq ($(MODE),RELEASE)
     FFLAGS = -cpp -O3 -mkl -assume bscc
     CHECKFFLAGS = -cpp -O0 -mkl -assume bscc
     
   else ifeq ($(MODE),VTUNE)
     FFLAGS = -g -debug -cpp -O3 -mkl -assume bscc -fopenmp
     CHECKFFLAGS = -g -debug -cpp -O0 -mkl -assume bscc -fopenmp

   endif
   LIBS += -L$(NETCDF_LIB) -lnetcdff $(PROBLEM_LIB)

endif

ifeq ($(COMM),PARALLEL)
  FFLAGS += -fopenmp
endif


$(info )
$(info )
# ------------- LIBRARIES ------------------------
INCLUDE= -I$(INC_DIR) -I$(NETCDF_INC) $(F90_INCLUDE)

# ------------- Main rules -----------------------
.PHONY: all	
.DEFAULT_GOAL := all

all: 	printHeader createFolders def_navier_stokes_macro ProblemFile $(EXECUTABLE)

runcyl: 	
	(cd ../TestCases/Cases/Cylinder && ../../../Solver/bin/$(EXECUTABLE) ./SETUP/Cylinder.HiOCase)		

debugcyl:
	(cd ../TestCases/Cases/FlowOverCircle && gdb ./$(EXECUTABLE))		

runvortex: 	
	(cd ../TestCases/Cases/TaylorVortex && ./$(EXECUTABLE) ./SETUP/TaylorVortex.HiOCase -check)		

runchan: 	
	(cd ../TestCases/Cases/ChannelFlow && ./$(EXECUTABLE) ./SETUP/Channel.HiOCase -check)		

runfree: 	
	(cd ../TestCases/Cases/FreeStream && ./$(EXECUTABLE) ./SETUP/FreeStream.HiOCase -check)		

Euler:  printHeader createFolders def_euler_macro ProblemFile $(EXECUTABLE)

def_euler_macro:
	$(eval EXECUTABLE:=HORSES2D.Euler)

def_navier_stokes_macro:
	$(eval NS_COMP:=-DNAVIER_STOKES)

# ---------- Executable rule -------------------
$(EXECUTABLE): $(OBJS) 
	$(FC) $(FFLAGS) $(NS_COMP) $(GLOB_MACROS) -o $(BIN_DIR)/$(EXECUTABLE) $(INCLUDE) \
	$(OBJS_EXEC) $(BUILD_DIR)/Checks.o $(LIBS)

# ---------- Library rule
ProblemFile: src/IO/ProblemFile.f90
	$(FC) -fPIC -shared $(FFLAGS) $(NS_COMP) $(GLOB_MACROS) $(INCLUDE) $(MOD_FLAG)$(INC_DIR) -c src/IO/ProblemFile.f90 -o $(BUILD_DIR)/$@.o
	$(FC) $(DYNLIB_FLAG) $(FFLAGS) $(NS_COMP) $(GLOB_MACROS) $(INCLUDE) $(MOD_FLAG)$(INC_DIR) ./build/ProblemFile.o

# ---------- Default compilation rule -----------
$(OBJS): % : $(or $(wildcard $(SRC_DIR)/**/%.f90), $(wildcard $(SRC_DIR)/%.f90))
	$(FC) $(FFLAGS) $(NS_COMP) $(GLOB_MACROS) $(INCLUDE) $(MOD_FLAG)$(INC_DIR) -c $(or $(wildcard $(SRC_DIR)/**/$@.f90),$(wildcard $(SRC_DIR)/$@.f90)) -o $(BUILD_DIR)/$@.o

# ---------- Rule for the Checks subroutine ----------
Checks: ./src/Generic/Checks.f90
	$(FC) $(CHECKFFLAGS) $(NS_COMP) $(GLOB_MACROS) $(INCLUDE) $(MOD_FLAG)$(INC_DIR) -c ./src/Generic/Checks.f90 -o $(BUILD_DIR)/Checks.o


clean:  printHeader
	$(info Removing installation:)
	rm -f -v -r $(BUILD_DIR)/*
	rm -f -v -r $(INC_DIR)/*
	rm -f -v -r $(LIB_DIR)/*

allclean: clean
	@rm -f -v -r $(BIN_DIR)/*

createFolders:
	@mkdir -p $(BIN_DIR) $(INC_DIR) $(LIB_DIR) $(BUILD_DIR)

printHeader:
	@echo '       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// '
	@echo '      |                                  _    _    ____    _____     _____   ______    _____                                   |'
	@echo '      |                                 | |  | |  / __ \  |  __ \   / ____| |  ____|  / ____|                                  |'
	@echo '      |                                 | |__| | | |  | | | |__) | | (___   | |__    | (___                                    |'
	@echo '      |                                 |  __  | | |  | | |  _  /   \___ \  |  __|    \___ \                                   |'
	@echo '      |                                 | |  | | | |__| | | | \ \   ____) | | |____   ____) |  				       |'
	@echo '      |                                 |_|  |_|  \____/  |_|  \_\ |_____/  |______| |_____/                                   |'
	@echo '      |                                                                                                                        |'
	@echo '      |                                  +m:y:-`                                       `-:y:m+                                 |'
	@echo '      |                          .:::odmMMMNMms:-                                   -:smMNMMMmdo:::.                           |'
	@echo '      |                          `+mMMNMMMMNNNd+:                                   :+dNNNMMMMNMMm+`                           |'
	@echo '      |                        /shdNMNymMMMMshMmdNmd:                           :dmNdmMhsMMMMmyNMNdhs/                         |'
	@echo '      |                       ooosMMMNMMMNmNhNMMMNmhd-                         -dhmNMMMNhNmNMMMNMMMsooo                        |'
	@echo '      |                      `-/+yhMMmMMMhhMhss+/+o+.                           .+o+/+sshMhhMMMmMMhy+/-`                       |'
	@echo '      |                          yMNMhMMMMNMs                                           sMNMMMMhMNMy                           |'
	@echo '      |                        :ydmMMhMMMMy+M/                                         /M+yMMMMhMMmdy:                         |'
	@echo '      |                       .+s+omMhMMMMMy+Ms`                                     `sM+yMMMMMhMmo+s+.                        |'
	@echo '      |                         +ymNMNmMMMMMmdMm+`                                 `+mMdmMMMMMmNMNmy+                          |'
	@echo '      |           :-           //oymMMmMMMMMMMmNMd-       /yys:       :syy/       -dMNmMMMMMMMmMMmyo//           -:            |'
	@echo '      |          - -y.          ./oymMMMMMMMMhddhNMNs..+yNMMMMm-     -mMMMMNy+..sNMNhddhMMMMMMMMmyo/.          .y- -           |'
	@echo '      |          +o.-N.       `-:+yhdNMMMMMsMMMMMMMMMMmMMMmd/msd     dsm/dmMMMmMMMMMMMMMMsMMMMMNdhy+:-`       .N-.o+           |'
	@echo '      |         : :NsMo           oyo+-hMMMMMMMMMMMMMM+NMd+`.+mM/   /Mm+.`+dMN+MMMMMMMMMMMMMMh-+oyo           oMsN: :          |'
	@echo '      |         /y.yMM:          `.   :hoMMMMMMMMMMMMMMmMNMMMMmsN` `NsmMMMMNMmMMMMMMMMMMMMMMoh:   .`          :MMy.y/          |'
	@echo '      |          oNMMh               oh+MMMMMMMMMMMMMMMMNdNNyNMdms smdMNyNNdNMMMMMMMMMMMMMMMM+ho               hMMNo           |'
	@echo '      |          oMMM-             /dyhMMMMMMMMMMMMMMMMNs/. `mmsMm mMsmm` ./sNMMMMMMMMMMMMMMMMhyd/             -MMMo           |'
	@echo '      |          dNym           -oNMMMMMMMMMMMMMMMMNoo/.    odh-mN Nm-hdo    ./ooNMMMMMMMMMMMMMMMMNo-           myNd           |'
	@echo '      |         `Mdsh        `+mMMMMMMMMMMMMMMMh+Mh-       -NM/ `. .` /MN-       -hM+hMMMMMMMMMMMMMMMm+`        hsdM`          |'
	@echo '      |         `NN:N.      oNMMMMMMNoMMMMMMMhoyh:       .oNMd`       `dMNo.       :hyohMMMMMMMoNMMMMMMNo      .N:NN`          |'
	@echo '      |          /Mm/do`  .mMMMMMMMMMMMMMMMmhyo.         NNM+           +MNN         .oyhmMMMMMMMMMMMMMMMm.  `od/mM/           |'
	@echo '      |           .smNMMdymMMMMMMMMNmMMMMMhy:            sy-             -ys            :yhMMMMMmNMMMMMMMMmydMMNms.            |'
	@echo '      |             -NhMMdMhmMMMMMMMyMMMhhMM+                                           +MMhhMMMyMMMMMMMmhMdMMhN-              |'
	@echo '      |             /yhNd-hNMMMMMMMMmymhmMMd                                             dMMmhmymMMMMMMMMNh-dNhy/              |'
	@echo '      |             :-sm/-`dMMMMMMMMhh/MMMm.                                             .mMMM/hhMMMMMMMMd`-/ms-:              |'
	@echo '      |               `yy  `oNMMMMMMMmsMMm.                                               .mMMsmMMMMMMMNo`  yy`                |'
	@echo '      |                 y`   `omMMMMMMsMy`    -/:/+++o.                       .o+++/:/-    `yMsMMMMMMmo`   `y                  |'
	@echo '      |                         sMMhmoNNhysyhmMMmNMMho`                       `ohMMNmMMmhysyhNNomhMMs                          |'
	@echo '      |                         oMMN:yMmyo/-.``                                       ``.-/oymMy:NMMo                          |'
	@echo '      |                        oNMN.                                                             .NMNo                         |'
	@echo '      |                        -sNmo.             High Order Spectral Element Solver            .omNs-                         |'
	@echo '      |                           /dyo.                                                       .oyd/                            |'
	@echo '      |                            `oMNd-      2D Compressible Navier-Stokes equations      -dNMo`                             |'
	@echo '      |                              :dMNy`                                               `yNMd:                               |'
	@echo '      |                               -ymms                                               smmy-                                |'
	@echo '      |                                                                                                                        |'
	@echo '       \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ '
	@echo ' '
	@echo ' '
	@echo ' '
	@echo ' *** Files: '
	@for i in $(OBJS_F90); do \
		echo "        > $$i ";\
	done
	@echo ' '
	@echo ' '
	@echo ' '
	@echo " *** Compiling with $(COMPILER) in $(MODE) mode"
	@echo ' '
	@echo ' '
	@echo ' '

