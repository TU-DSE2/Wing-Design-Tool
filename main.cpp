#include <iostream>

#include <armadillo>

#include "CFDSolver/CFDSolver.h"

using namespace std;
using namespace arma;

int main() {
    CFDSolver cfdsolver;
    cfdsolver.addVortex(cfdsolver.Nx*3/16, cfdsolver.Ny*1/2, (double)cfdsolver.Nx/8, 1);
    cfdsolver.run(1);
    cout << cfdsolver.w;
}
