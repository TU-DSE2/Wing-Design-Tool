#ifndef GRIDBDRY_H_INCLUDED
#define GRIDBDRY_H_INCLUDED

#include <armadillo>

using namespace arma;

/*
Created:        4 december 2013     by Jeroen Barnhoorn
Last update:    4 december 2013     by Jeroen Barnhoorn

STATUS: WIP
*/
class GridBdry {
    public:
        mat Gamma;
        mat n;
        mat OmegaI;

        GridBdry(mat, mat, mat);

        void setGamma(mat);
        void setN(mat);
        void setOmegaI(mat);
        int npoints();
        rowvec interp_normal(cube u);     //TODO
};

#endif // GRIDBDRY_H_INCLUDED
