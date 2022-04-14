# Computational-Magnetohydrodynamics
This directory contains all the code required to perform 1-D and 2-D simulations of magnetohydrodynamics system. Final report `MHD_report.pdf` submitted has also been added.

* **main.cpp** - This file contains the initial settings required to run simulations.
* **array.cpp** - This file contains the functions defining the data structure used in the MHD code.
* **array.H** - This file contains the function headers for `array.cpp`.
* **system_MHD.cpp** - This file contains the functions for computing time step, setting boundary condition, performing iteration and setting initial data for different MHD test cases
* **system_MHD.H** - This file contains the function headers for `system_MHD.cpp`.
* **eos_MHD.cpp** - This file contains the functions for MHD specific flux functions and variable conversion.
* **eos_MHD.H** - This file contains the function headers for `eos_MHD.cpp`.
* **numerical_method_MHD.cpp** - This file contains the functions for SLIC numerical solver.
* **numerical_method_MHD.H** - This file contains the function headers for `numerical_method_MHD.cpp`.
* **convergence.cpp** - This file contains the code for performing the convergence analysis on the simulation results. 
* **RiemannExactSol.cpp** - This file contains the solver for computing the exact results of 1-D Euler.


