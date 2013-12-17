#ifndef PROFILE_H
#define PROFILE_H

#include <armadillo>
#include "GridBdry.h"
#include "ExactBdry.h"

using namespace arma;

/*
Created:        13 december 2013    by Jeroen Barnhoorn
Last updated:   13 december 2013    by Jeroen Barnhoorn

STATUS: DONE
*/

class Profile {
    public:
        Profile();
        GridBdry getGridBdry();
        ExactBdry getExactBdry();
    protected:
        mat Gamma_I;
        mat n_grid;
        mat interior;
        mat ex_loc;
        mat ex_n;
        vec ex_ds;
};

#endif // PROFILE_H
