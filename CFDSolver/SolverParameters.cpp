#include "include/SolverParameters.h"

/*
  Created:         3 december 2013     by Jeroen Barnhoorn
  Last update:    11 december 2013     by Jeroen Barnhoorn

  STATUS: DONE

  Sets solver parameters
*/

SolverParameters::SolverParameters() {
    boundaries   = true;                    //Use boundaries
    viscous      = true;                    //Use viscous flow
    viscosity    = 1e-2;                    //Viscosity value
    mean_flow[0] = 1;                       //Far flow speed (x-direction)
    mean_flow[1] = -0;                    //Far flow speed (y-direction)
    deltat       = -1;

    grid_size[0]                = 32;       //Number of grid points (x-direction)
    grid_size[1]                = 32;       //Number of grid points (y-direction)
    n_advection_steps           = 1;
    max_iter_streamfn           = 2;
    min_residual_streamfn       = 1e-4;
    max_iter_potential          = 3;
    min_residual_potential      = 1e-3;
    max_iter_potential_setup    = 50;
    min_residual_potential_setup = 1e-12;
    max_iter_GMRes              = 5;
    min_residual_GMRes          = 1e-3;

    matplotlib_output           = -1;

    use_GPU     = true;                     //Use GPU for expensive calculations
}
