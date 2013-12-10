#ifndef CFDSOLVER_H
#define CFDSOLVER_H

#include <armadillo>

#include "include/SolverParameters.h"
#include "include/Rod.h"
#include "include/GridBdry.h"
#include "include/ExactBdry.h"
#include "include/PoissonSolver.h"
#include "include/DiffOps.h"
#include "include/Interpolate.h"
#include "include/viscous_wall_vorticity_flux.h"

using namespace std;
using namespace arma;

class CFDSolver {
    public:
        CFDSolver();

        GridBdry getGridBdry();
        void run(int);
        void add_vortex(int, int, double, double);
        float** getWArray();

        int Nx;
        int Ny;
        float nu;
        cube u;
        mat w;
        mat ddw;
        mat psi;
        mat phi;
        cube dphi;
        mat g_force;

        SolverParameters solverparameters;
        GridBdry gridbdry;
        ExactBdry exactbdry;
        Rod Profile;
        vec x;
        vec y;
        vec g;
        rowvec garg;
        rowvec b;
        mat Lmat;
        double deltat;
        int iter;
        double t;
    private:
        rowvec L_matvec(rowvec);
        mat L_matrix();
        int iterations;
        int total_iterations;
        cube abs_u;
};

#endif // CFDSOLVER_H
