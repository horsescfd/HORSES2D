COMPILER ?=gfortran
INC_DIR=./include
BUILD_DIR=./build
BIN_DIR=./bin
SRC_DIR=./src
EXECUTABLE = HiODG2D
MODE ?=DEBUG
$(info $(MODE))
BLAS_LIB = 
LAPACK_LIB = 
LIBS=

# --------------- MACROS ---------------------------------------------------

GLOB_MACROS ?= -D_USE_LAPACK

# --------------- DEPENDENCIES ---------------------------------------------------


Physics: PhysicsNS Headers
PhysicsNS: Setup SMConstants Headers
DGSEM_class: Physics SMConstants NodesAndWeightsClass Setup Mesh1D MeshFile DGSpatialDiscretizationMethods DGTimeIntegrator Storage DGBoundaryConditions
NodesAndWeightsClass: SMConstants InterpolationAndDerivatives LegendreAlgorithms MatrixOperations
InterpolationAndDerivatives: SMConstants
LegendreAlgorithms: SMConstants
Mesh1D: FaceClass NodeClass Element1DClass MeshFile Setup InitialConditions Storage
MeshFile: Setup SMConstants
Element1DClass: NodeClass SMConstants Physics NodesAndWeightsClass Storage
NodeClass: SMConstants
FaceClass: Element1DClass NodeClass SMConstants
Checks: DGSEM_class Physics SMConstants NodesAndWeightsClass Setup Mesh1D MeshFile DGSpatialDiscretizationMethods Storage DGBoundaryConditions
Headers: SMConstants
NetCDF_Interface: SMConstants
ProblemFile: SMConstants Setup
InitialConditions: SMConstants Physics
Setup: SMConstants
DGSpatialDiscretizationMethods: SMConstants Physics Mesh1D DGSecondOrderMethods MatrixOperations DGFirstOrderMethods
DGSecondOrderMethods: SMConstants Physics Setup FaceClass Element1DClass MatrixOperations NodesAndWeightsClass
DGFirstOrderMethods: SMConstants Physics Setup FaceClass Element1DClass MatrixOperations NodesAndWeightsClass
DGBoundaryConditions: SMConstants FaceClass Setup Physics
LAPackInterface: SMConstants
FileWriting: SMConstants Physics Mesh1D Setup
DGTimeIntegrator: SMConstants Mesh1D DGSpatialDiscretizationMethods Setup FileWriting
Storage: SMConstants
main: SMConstants NodesAndWeightsClass Physics MeshFile DGSEM_class


# --------------- DO NOT EDIT ----------------------------------------------------

OBJS_F90 = $(notdir $(wildcard $(SRC_DIR)/**/*.f90)) $(notdir $(wildcard $(SRC_DIR)/*.f90))
OBJS= $(patsubst %.f90,%,$(OBJS_F90))
OBJS_EXEC = $(foreach obj,$(OBJS),$(BUILD_DIR)/$(obj).o) 


include make.inc

# -------------- CHOOSE COMPILER ------------------
ifeq ($(COMPILER),gfortran)
   MOD_FLAG :=-J
   FC=gfortran
   ifeq ($(MODE),DEBUG)
      $(info MODE = DEBUG)
      FFLAGS = -ffree-line-length-0 -cpp -O0 -g -fcheck=all -Wno-unused-variable -fbacktrace -fbounds-check -finit-real=snan -Wall -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries -fbackslash -fcheck=pointer
   else 
  	$(info MODE = RELEASE)
      FFLAGS = -cpp -ffree-line-length-0 -O3 -ftree-vectorize  -ftree-vectorizer-verbose=0 -fbackslash
   endif
   LIBS += -L$(BLAS_LIB) -L$(LAPACK_LIB) -llapack -lblas -L$(NETCDF_LIB) -lnetcdff

else ifeq ($(COMPILER),ifort)
   MOD_FLAG := -module 
   FC=ifort
   ifeq ($(MODE),DEBUG)
     $(info MODE = DEBUG)
     FFLAGS = -cpp -O0 -g -warn all -fpscomp logicals -mkl -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces -ftrace=full -assume ieee_fpe_flags -ftz -fp-model precise -fp-speculation=off -assume bscc
   else ifeq ($(MODE),RELEASE)
     $(info MODE = RELEASE)
     FFLAGS = -cpp -O3 -mkl -assume bscc
     
   endif
   LIBS += -L$(NETCDF_LIB) -lnetcdff

endif


# ------------- LIBRARIES ------------------------
INCLUDE= -I$(INC_DIR) -I$(NETCDF_INC)

# ------------- Main rules -----------------------
	

all: 	$(EXECUTABLE)
	
# ---------- Executable rule -------------------
$(EXECUTABLE): $(OBJS) 
	$(FC) $(FFLAGS) -o $(BIN_DIR)/$(EXECUTABLE) $(INCLUDE) \
	$(OBJS_EXEC) $(LIBS)\

# ---------- Default compilation rule -----------
$(OBJS): % : $(or $(wildcard $(SRC_DIR)/**/%.f90), $(wildcard $(SRC_DIR)/%.f90))
	$(FC) $(FFLAGS) $(GLOB_MACROS) $(INCLUDE) $(MOD_FLAG)$(INC_DIR) -c $(or $(wildcard $(SRC_DIR)/**/$@.f90),$(wildcard $(SRC_DIR)/$@.f90)) -o $(BUILD_DIR)/$@.o


clean:
	rm -f $(BUILD_DIR)/*
	rm -f $(INC_DIR)/*

allclean: clean
	@rm -r -f $(BIN_DIR)/*

