#ifndef ROD_H
#define ROD_H

#include <armadillo>
#include "GridBdry.h"
#include "ExactBdry.h"

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
        Rod(vec, vec);

        GridBdry getGridBdry();
        ExactBdry getExactBdry();
    protected:
    private:
        mat Gamma_I;
        mat n;
        mat interior;
        mat ex_loc;
        mat ex_n;
        vec ex_ds;
};

#endif // ROD_H
