#ifndef CFDSOLVER_H
#define CFDSOLVER_H

#include <armadillo>

#include "include/SolverParameters.h"
#include "include/Rod.h"
#include "include/Cylinder.h"
#include "include/CustomProfile.h"
#include "include/GridBdry.h"
#include "include/ExactBdry.h"
#include "include/PoissonSolver.h"
#include "include/DiffOps.h"
#include "include/Interpolate.h"
#include "include/ViscousWallVorticityFlux.h"

using namespace std;
using namespace arma;

class CFDSolver {
    public:
        CFDSolver();
        CFDSolver(int, int);

        GridBdry getGridBdry();
        void run(int);
        void addVortex(int, int, double, double);

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
        //CustomProfile profile;
        //Cylinder profile;
        Rod profile;
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
