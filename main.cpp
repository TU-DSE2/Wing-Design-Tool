#include <iostream>

#include <armadillo>
//#include <cmath>

#include "include/SolverParameters.h"
#include "include/DiffOps.h"
#include "include/PoissonSolver.h"
#include "include/Interpolate.h"
//#include "include/Geometries.h"
#include "include/Rod.h"

//using namespace std;
using namespace arma;

bool verbose = true;

int main() {
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

    Rod Profile(x, y);
    GridBdry gridbdry = Profile.getGridBdry();
    ExactBdry exactbdry = Profile.getExactBdry();

    return 0;
}
