#include "../include/GridBdry.h"

/*
Created:         4 december 2013     by Jeroen Barnhoorn
Last update:    10 december 2013    by Jeroen Barnhoorn

STATUS: WIP
*/

GridBdry::GridBdry(mat loc, mat interior, mat nn) {
    Gamma = loc;
    OmegaI = interior;
    n_grid = nn;
}

void GridBdry::setPoints(mat loc, mat interior, mat nn) {
    Gamma = loc;
    n_grid = nn;
    OmegaI = interior;
}

void GridBdry::setGamma(mat loc) {
    Gamma = loc;
}

void GridBdry::setN(mat n) {
    n_grid = n;
}

void GridBdry::setOmegaI(mat interior) {
    OmegaI = interior;
}

int GridBdry::npoints() {
    return Gamma.n_rows;
}

rowvec GridBdry::interp_normal(cube u) {      //TODO
    rowvec temp = zeros<rowvec>(Gamma.n_rows);
    for(unsigned int i = 0; i < Gamma.n_rows; i++) {
        temp(i) = u(Gamma(i, 0), Gamma(i, 1), 0)*n_grid(i, 0) + u(Gamma(i, 0), Gamma(i, 1), 1)*n_grid(i, 1);
    }
    return temp;
}
