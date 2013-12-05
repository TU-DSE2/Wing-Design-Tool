#include <iostream>

#include <armadillo>
//#include <cmath>

#include "include/SolverParameters.h"
#include "include/DiffOps.h"
#include "include/PoissonSolver.h"
#include "include/Interpolate.h"
//#include "include/Geometries.h"
#include "include/Rod.h"
#include "include/Vortex.h"
#include "include/Func.h"

//using namespace std;
using namespace arma;

SolverParameters SolverParameters;

int Nx = SolverParameters.grid_size[0];
int Ny = SolverParameters.grid_size[1];
vec x = linspace<vec>(0, Nx-1, Nx);
vec y = linspace<vec>(0, Ny-1, Ny);
float nu = SolverParameters.viscosity;
int iter = 0;
double t = 0;
float source = 0.1;

cube u = zeros<cube>(Nx, Ny, 2);
mat w = zeros<mat>(Nx, Ny);
mat ddw = zeros<mat>(Nx, Ny);
mat psi = zeros<mat>(Nx, Ny);
mat phi = zeros<mat>(Nx, Ny);
cube dphi = zeros<cube>(Nx, Ny, 2);
mat g_force = zeros<mat>(Nx, Ny);

rowvec L_matvec(rowvec garg, GridBdry gridbdry) {
    for(int i = 0; i < garg.n_cols; i++) {
        g_force(gridbdry.Gamma(i, 0), gridbdry.Gamma(i, 1)) = garg(i);
    }
    phi = pSolve(g_force);
    dphi = grad_bdry(phi, gridbdry);
    return gridbdry.interp_normal(dphi);
}

mat L_matrix(GridBdry gridbdry) {
    int M = gridbdry.npoints();
    mat Lmat = zeros<mat>(M, M);
    rowvec garg = zeros<rowvec>(M);
    for(unsigned int i = 0; i < M; i++) {
        garg = zeros<rowvec>(M);
        garg(i) = 1;
        L_matvec(garg, gridbdry);
        Lmat.col(i) = trans(L_matvec(garg, gridbdry));
    }
    return Lmat;
}

int main() {
    Rod Profile(x, y);
    GridBdry gridbdry = Profile.getGridBdry();
    ExactBdry exactbdry = Profile.getExactBdry();
    vec g = zeros<vec>(gridbdry.npoints());

    w = add_vortex(w, x, y, Nx*3/16, Ny*1/2, Nx/8, 1);

    cout << sum(sum(L_matrix(gridbdry)));

    return 0;
}
