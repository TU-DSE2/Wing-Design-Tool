#include "include/SolverParameters.h"

/*
Created:        3 december 2013     by Jeroen Barnhoorn
Last update:    4 december 2013     by Jeroen Barnhoorn

STATUS: DONE
*/
SolverParameters::SolverParameters() {
    boundaries   = true;
    viscous      = true;
    viscosity    = 1e-2;
    mean_flow[0] = 1;
    mean_flow[1] = -0.5;
    deltat       = -1;

    grid_size[0]            = 9;
    grid_size[1]            = 9;
    n_advection_steps       = 1;
    max_iter_streamfn       = 2;
    min_residual_streamfn   = 1e-4;
    max_iter_potential      = 3;
    min_residual_potential  = 1e-3;
    max_iter_potential_setup = 50;
    min_residual_potential_setup = 1e-12;
    max_iter_GMRes          = 5;
    min_residual_GMRes      = 1e-3;

    matplotlib_output       = -1;
}
