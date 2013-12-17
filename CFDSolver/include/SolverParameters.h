#ifndef SOLVERPARAMETERS_H
#define SOLVERPARAMETERS_H

class SolverParameters {
    public:
        SolverParameters();

        bool boundaries;
        bool viscous;
        float viscosity;
        float mean_flow [2];
        float deltat;

        int grid_size [2];
        int n_advection_steps;
        int max_iter_streamfn;
        float min_residual_streamfn;
        int max_iter_potential;
        float min_residual_potential;
        int max_iter_potential_setup;
        float min_residual_potential_setup;
        int max_iter_GMRes;
        float min_residual_GMRes;

        int matplotlib_output;

        bool use_cuda;
};

#endif // SOLVERPARAMETERS_H
