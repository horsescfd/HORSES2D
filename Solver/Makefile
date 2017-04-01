$(info High-Order Discontinuous Galerkin 2D Solver Compilation. Madrid TU.)
$(info )
$(info )


COMPILER ?=gfortran
INC_DIR=./include
BUILD_DIR=./build
BIN_DIR=./bin
SRC_DIR=./src
EXECUTABLE = HiODG2D
MODE ?=DEBUG
BLAS_LIB = 
LAPACK_LIB = 
LIBS=
GLOB_MACROS := -D_DIMENSIONLESS_TAU

# --------------- MACROS ---------------------------------------------------

#GLOB_MACROS +=-D_USE_LAPACK

# --------------- DEPENDENCIES ---------------------------------------------------


Physics: PhysicsNS Headers
PhysicsNS: Setup SMConstants Headers MatrixOperations
VariableConversion: PhysicsNS
RiemannSolvers: PhysicsNS SMConstants InviscidFluxes
InviscidFluxes: PhysicsNS SMConstants
ViscousFluxes: PhysicsNS SMConstants
DGSEM_class: Physics SMConstants NodesAndWeightsClass Setup QuadMesh MeshFile DGSpatialDiscretizationMethods DGTimeIntegrator Storage DGBoundaryConditions Tecplot
NodesAndWeightsClass: SMConstants InterpolationAndDerivatives LegendreAlgorithms MatrixOperations Physics
InterpolationAndDerivatives: SMConstants Utilities
LegendreAlgorithms: SMConstants
QuadMesh: NodeClass QuadElement MeshFile Setup InitialConditions Storage DGBoundaryConditions Utilities
MeshFile: Setup SMConstants NetCDF_Interface Headers
QuadElement: NodeClass SMConstants Physics NodesAndWeightsClass Storage InterpolationAndDerivatives MatrixOperations MetricMatrix
NodeClass: SMConstants Physics
MetricMatrix: SMConstants
Checks: DGSEM_class Physics SMConstants NodesAndWeightsClass Setup QuadMesh MeshFile DGSpatialDiscretizationMethods Storage DGBoundaryConditions
Headers: SMConstants
NetCDF_Interface: SMConstants
ProblemFile: SMConstants Setup
InitialConditions: SMConstants Physics
Setup: SMConstants Paramfile
DGWeakIntegrals: SMConstants Physics QuadElement MatrixOperations
DGSpatialDiscretizationMethods: SMConstants Physics QuadMesh DGViscousMethods MatrixOperations DGInviscidMethods DGWeakIntegrals
DGViscousMethods: SMConstants Physics Setup QuadElement MatrixOperations NodesAndWeightsClass QuadMesh Headers
DGInviscidMethods: SMConstants Physics Setup QuadElement MatrixOperations NodesAndWeightsClass
DGInviscid_StandardDG: DGInviscidMethods
DGViscous_BR1: DGViscousMethods
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
DGTimeIntegrator: SMConstants QuadMesh DGSpatialDiscretizationMethods Setup FileWriting Tecplot Monitors
Monitors: SMConstants Setup Paramfile Physics QuadMesh
Storage: SMConstants
main: SMConstants NodesAndWeightsClass Physics MeshFile DGSEM_class


# --------------- DO NOT EDIT ----------------------------------------------------

F90_REMOVE=
F90_INCLUDE=-I$(SRC_DIR)/Mesh/ -I$(SRC_DIR)/Generic/

OBJS_F90:= $(notdir $(wildcard $(SRC_DIR)/**/*.f90)) $(notdir $(wildcard $(SRC_DIR)/*.f90))
OBJS_F90:= $(filter-out $(F90_REMOVE),$(OBJS_F90))
$(info )
$(info ***Files: $(OBJS_F90) )
$(info )
OBJS= $(patsubst %.f90,%,$(OBJS_F90))
OBJS_EXEC = $(foreach obj,$(OBJS),$(BUILD_DIR)/$(obj).o) 


include make.inc

# -------------- CHOOSE COMPILER ------------------
ifeq ($(COMPILER),gfortran)
   MOD_FLAG :=-J
   FC=gfortran
   ifeq ($(MODE),DEBUG)
      $(info Compiling with gfortran in DEBUG mode.)
#      FFLAGS = -ffree-line-length-0 -cpp -O0 -g -fcheck=all -Wno-unused-variable -fbacktrace -fbounds-check -finit-real=snan -Wall -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries -fbackslash -fcheck=pointer
       FFLAGS = -cpp -O0 -g -fcheck=all -ffree-line-length-0 -fbacktrace -fbounds-check -Wall -fcheck-array-temporaries -fbackslash -ffpe-trap=invalid,zero,overflow -static-libgfortran
   else 
  	$(info Compiling with gfortran in RELEASE mode.)
      FFLAGS = -cpp -ffree-line-length-0 -O3 -ftree-vectorize  -ftree-vectorizer-verbose=0 -fbackslash
   endif
   LIBS += -L$(BLAS_LIB) -L$(LAPACK_LIB) -llapack -lblas -L$(NETCDF_LIB) -lnetcdff

else ifeq ($(COMPILER),ifort)
   MOD_FLAG := -module 
   FC=ifort
   ifeq ($(MODE),DEBUG)
     $(info Compiling with intel fortran in DEBUG mode.)
     FFLAGS = -cpp -O0 -g -warn all -fpscomp logicals -mkl -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces -ftrace=full -assume ieee_fpe_flags -ftz -fp-model precise -fp-speculation=off -assume bscc
   else ifeq ($(MODE),RELEASE)
     $(info Compiling with intel fortran in RELEASE mode.)
     FFLAGS = -cpp -O3 -mkl -assume bscc
     
   else ifeq ($(MODE),VTUNE)
     $(info Compiling with intel fortran in VTUNE mode.)
     FFLAGS = -g -debug -cpp -O3 -mkl -assume bscc

   endif
   LIBS += -L$(NETCDF_LIB) -lnetcdff

endif


$(info )
$(info )
# ------------- LIBRARIES ------------------------
INCLUDE= -I$(INC_DIR) -I$(NETCDF_INC) $(F90_INCLUDE)

# ------------- Main rules -----------------------
	

all: 	def_navier_stokes_macro $(EXECUTABLE)

runcyl: 	
	(cd ../TestCases/Benchmarktests/FlowOverCircle && ./$(EXECUTABLE) ./SETUP/Cylinder.HiOCase -check)		

debugcyl:
	(cd ../TestCases/Benchmarktests/FlowOverCircle && gdb ./$(EXECUTABLE))		

runvortex: 	
	(cd ../TestCases/Benchmarktests/TaylorVortex && ./$(EXECUTABLE) ./SETUP/TaylorVortex.HiOCase -check)		

runchan: 	
	(cd ../TestCases/Benchmarktests/ChannelFlow && ./$(EXECUTABLE) ./SETUP/Channel.HiOCase -check)		

runfree: 	
	(cd ../TestCases/Benchmarktests/FreeStream && ./$(EXECUTABLE) ./SETUP/FreeStream.HiOCase -check)		

Euler:  def_euler_macro $(EXECUTABLE)

def_euler_macro:
	$(eval EXECUTABLE:=HiODG2D.Euler)

def_navier_stokes_macro:
	$(eval NS_COMP:=-DNAVIER_STOKES)

# ---------- Executable rule -------------------
$(EXECUTABLE): $(OBJS) 
	$(FC) $(FFLAGS) $(NS_COMP) $(GLOB_MACROS) -o $(BIN_DIR)/$(EXECUTABLE) $(INCLUDE) \
	$(OBJS_EXEC) $(LIBS)\

# ---------- Default compilation rule -----------
$(OBJS): % : $(or $(wildcard $(SRC_DIR)/**/%.f90), $(wildcard $(SRC_DIR)/%.f90))
	$(FC) $(FFLAGS) $(NS_COMP) $(GLOB_MACROS) $(INCLUDE) $(MOD_FLAG)$(INC_DIR) -c $(or $(wildcard $(SRC_DIR)/**/$@.f90),$(wildcard $(SRC_DIR)/$@.f90)) -o $(BUILD_DIR)/$@.o


clean:
	$(info Removing installation:)
	rm -f $(BUILD_DIR)/*
	rm -f $(INC_DIR)/*

allclean: clean
	@rm -r -f $(BIN_DIR)/*

