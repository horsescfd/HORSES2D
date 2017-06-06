# HORSES: a High Order Spectral Element Solver

![image](https://dl.dropboxusercontent.com/s/kj8zqel72zyolgv/Logo1.png?dl=0)


***This is a two dimensional discontinuous Galerkin spectral element method solver for Navier-Stokes equations.***

This code has been developed in the Madrid Technical University.

## **Dependencies**

### Operating system

This code has been tested and its currently supported in:

  * MacOSX (OSX Sierra)
  * Linux (Ubuntu 14.04)
  
It still has not been tested in Windows OS.

### FORTRAN Compiler

This code has been tested with:

  * GNU Fortran version 6.3.0: https://github.com/gcc-mirror/gcc.git
  * Intel Fortran 2017.

Earlier versions may not be supported, since the code uses fortran submodules, which is a recent fortran 2008 feature.

### NetCDF Fortran Libraries

Network Common Data (NetCDF) libraries are required to run the code. This libraries are used for high performance data files processing, such as mesh and results files. Installing NetCDF libraries can be tricky for the first time. You can find all the information in their website:

https://www.unidata.ucar.edu/downloads/netcdf/index.jsp

If you attempt a manual installation, recall that you will need this prerequisites to install NetCDF-C libraries.

  * HDF5 1.8.9 or later (for netCDF-4 support)
  * zlib 1.2.5 or later (for netCDF-4 compression)
  * curl 7.18.0 or later (for DAP remote access client support)

### BLAS and LAPACK

BLAS (Basic Linear Algebra Solvers) and LAPack (Linear Algebra Package) are optional, and they are just used for matricial operations. This is controlled by the flag

```
-D_USE_LAPACK
```

in the code Makefile.

### Postprocessing tools

This code supports either **Tecplot** and **Paraview**. Using the former or the latter is selected by a runtime variable in each job case file.

## Description:

This code is divided in three parts:

1. Common binaries, which are stored in the ./Solver/bin folder. Two binaries are generated:
  * **HORSES2D.NS**: The **Navier-Stokes** equations solver.
  * **HORSES2D.Euler**: The **Euler** equations solver.

2. A shared library, named "**libproblemfile** which contain several case-dependant subroutines. A default **libproblemfile** library is compiled in **./Solver/lib**, which can be replaced in runtime by case specific libraries. During runtime, the binaries will link the library in this order:
  1. In  **./SETUP/libproblemfile**: Case specific libproblemfile.
  2. In **$(HORSES_PATH)/Solver/lib**: Default libproblemfile
  3. In **$(LD_LIBRARY_PATH/DYLD_LIBRARY_PATH)**: Additional path.
  
3. A case file which controls the simulation parameters. An example can be found [here](Utils/CaseFile/DefaultCaseFile.HiOCase). Flow and boundary conditions, solver parameters, and simulation configurations are selected in this file.

## Compilation

The code compilation uses GNU make (https://www.gnu.org/software/make/). The library dependencies (such as NetCDF, BLAS, and LAPACK) are speficied in a **make.inc** file. Once this file is placed in **./Solver**, the code is ready to be compiled:


```
	>> ./configure.sh  (or sh ./configure.sh if no executable permissions are given to the configure.sh script)
 	>> cd ./Solver
	>> make allclean
	>> make COMPILER=gfortran/ifort MODE=DEBUG/RELEASE COMM=SERIAL/PARALLEL 
	>> make COMPILER=gfortran/ifort MODE=DEBUG/RELEASE COMM=SERIAL/PARALLEL Euler
```


This creates the binaries in the **./Solver/bin** folder, and the default **libproblemfile** library in the **./Solver/lib** folder.

Next, to compile a test case specific libproblemfile library, a Makefile.template can be found in each test case **./SETUP** folder, for instance:

```
	>> cd $(TEST_CASE)/SETUP
	>> make -f Makefile.template MODE=RELEASE/DEBUG COMPILER=ifort/gfortran COMM=SERIAL/PARALLEL
```

## Execution

The binary will look for the case file in the first command line argument, that is:

```
	>> HORSES2D.NS Casefile.HiOCase
```
  
