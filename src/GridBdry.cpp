#include "../include/GridBdry.h"

/*
Created:        4 december 2013     by Jeroen Barnhoorn
Last update:    5 december 2013     by Jeroen Barnhoorn

STATUS: WIP
*/

GridBdry::GridBdry(mat loc, mat interior, mat nn) {
    Gamma = loc;
    n = nn;
    OmegaI = interior;
}

void GridBdry::setGamma(mat loc) {
    Gamma = loc;
}

void GridBdry::setN(mat n) {
    n = n;
}

void GridBdry::setOmegaI(mat interior) {
    OmegaI = interior;
}

int GridBdry::npoints() {
    return Gamma.n_rows;
}

rowvec GridBdry::interp_normal(cube u) {      //TODO
    rowvec temp = zeros<rowvec>(Gamma.n_rows);
    for(int i = 0; i < Gamma.n_rows; i++) {
        temp(i) = u(Gamma(i, 0), Gamma(i, 1), 0)*n(i, 0) + u(Gamma(i, 0), Gamma(i, 1), 1)*n(i, 1);
    }
    return temp;
}
