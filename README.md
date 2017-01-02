This is a two dimensional discontinuous Galerkin spectral element method solver for Navier-Stokes equations. This code has been developed in the Madrid Technical University.

Installation:
	
	1) Replace make.inc.example to make.inc with the paths of your system.
	2) Run make all. Several options are available:
		- Mode DEBUG/RELEASE
		- Compiler gfortran/ifort

	3) Run the three benchmark tests:
		- make runcyl. Flow around a circle
		- make runvortex. Taylor vortex problem.
		- make runchan. Channel.
