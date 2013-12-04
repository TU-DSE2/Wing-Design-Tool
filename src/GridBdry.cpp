#include "../include/GridBdry.h"

/*
Created:        4 december 2013     by Jeroen Barnhoorn
Last update:    4 december 2013     by Jeroen Barnhoorn

STATUS: WIP
*/
GridBdry::GridBdry() {
}

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

vec GridBdry::interp_normals(cube u) {      //TODO
    vec temp = zeros<vec>(Gamma.n_rows);
    return temp;
}
