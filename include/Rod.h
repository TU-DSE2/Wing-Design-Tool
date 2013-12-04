#ifndef ROD_H
#define ROD_H

#include <armadillo>
#include "GridBdry.h"

using namespace arma;

/*
Written by: Jeroen Barnhoorn,
Created: 4 december 2013

STATUS: WIP
*/
class Rod {
    public:
        vec x;
        vec y;
        mat Gamma_I;
        mat n;
        mat interior;

        Rod(vec, vec);

        GridBdry getGridBdry();
    protected:
    private:
};

#endif // ROD_H
