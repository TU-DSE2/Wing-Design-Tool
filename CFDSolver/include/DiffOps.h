#ifndef DIFFOPS_INCLUDED
#define DIFFOPS_INCLUDED

#include <armadillo>

#include "GridBdry.h"

using namespace arma;

/*
Created:         4 december 2013    by Steve Brust
Last updated:   10 december 2013    by Jeroen Barnhoorn

STATUS: DONE
*/
static mat laplace(mat w) {
    mat ddw = zeros<mat>(w.n_rows, w.n_cols);

    ddw.cols(1, ddw.n_cols - 2) = w.cols(2, w.n_cols - 1)
    		- 2*w.cols(1, w.n_cols - 2) + w.cols(0, w.n_cols - 3); //indices take
    //into account actual width and index start at zero. Middle cols

    ddw.rows(1, ddw.n_rows - 2) = ddw.rows(1, ddw.n_rows - 2)
    		+ w.rows(2, w.n_rows - 1) - 2*w.rows(1, w.n_rows - 2)
    		+ w.rows(0, w.n_rows - 3); //+= statement for middle rows

    ddw.col(0) = ddw.col(0) + 2*(w.col(1) - w.col(0)); //+= statement on first col ,
    //Farfield Boundaries

    ddw.col(ddw.n_cols - 1) = ddw.col(ddw.n_cols - 1) + 2*(w.col(w.n_cols - 2)
    		- w.col(w.n_cols - 1)); //+= statement last col

    ddw.row(0) = ddw.row(0) + 2*(w.row(1) - w.row(0));//+= statement on first row//GOOD

    ddw.row(ddw.n_rows - 1) = ddw.row(ddw.n_rows - 1)
    		+ 2*(w.row(w.n_rows - 2) - w.row(w.n_rows - 1));//+= last row //GOOD
    return ddw;
}

/*
Created:        18 november 2013    by Jeroen Barnhoorn
Last updated:   10 december 2013    by Jeroen Barnhoorn

STATUS: DONE
*/
static cube grad(mat phi) {
    cube u = zeros<cube>(phi.n_rows, phi.n_cols, 2);
    for(unsigned int i = 0; i < phi.n_rows; i++) {
        for(unsigned int j = 0; j < phi.n_cols; j++) {
            if(j > 0 && j < (phi.n_cols-1)) {
                u(i, j, 1) = 0.5*(phi(i, j+1) - phi(i, j-1));
            }
            if(i > 0 && i < (phi.n_rows-1)) {
                u(i, j, 0) = 0.5*(phi(i+1, j) - phi(i-1, j));
            }
            if(i == 0) {
                u(i, j, 0) = phi(i+1, j) - phi(i, j);
            }
            if(i == (phi.n_rows-1)) {
                u(i, j, 0) = phi(i, j) - phi(i-1, j);
            }
            if(j == 0) {
                u(i, j, 1) = phi(i, j+1) - phi(i, j);
            }
            if(j == phi.n_cols-1) {
                u(i, j, 1) = phi(i, j) - phi(i, j-1);
            }
        }
    }
    return u;
}

/*
Created:        18 november 2013    by Jeroen Barnhoorn
Last updated:   10 december 2013    by Jeroen Barnhoorn

Curl of vector => scalar

STATUS: DONE
*/
static mat curl(cube u) {
    mat w = zeros(u.n_rows, u.n_cols);
    for(unsigned int i = 0; i < u.n_rows; i++) {
        for(unsigned int j = 0; j < u.n_cols; j++) {
            if(i > 0 && i < (u.n_cols-1)) {
                w(i, j) = 0.5*(u(i+1, j, 1) - u(i-1, j, 1));
            }
            if(j > 0 && j < (u.n_rows-1)) {
                w(i, j) = w(i, j) - 0.5*(u(i, j+1, 0) - u(i, j-1, 0));
            }
            if(i == 0) {
                w(i, j) = w(i, j) + u(i+1, j, 1) - u(i, j, 1);
            }
            if(i == u.n_cols-1) {
                w(i, j) = w(i, j) + u(i, j, 1) - u(i-1, j, 1);
            }
            if(j == 0) {
                w(i, j) = w(i, j) - u(i, j+1, 0) + u(i, j, 0);
            }
            if(j == u.n_rows-1) {
                w(i, j) = w(i, j) - u(i, j, 0) + u(i, j-1, 0);
            }
        }
    }
    return w;
}

/*
Created:        19 november 2013    by Jeroen Barnhoorn
Last changed:   10 december 2013    by Jeroen Barnhoorn

Curl of scalar => vector

STATUS: DONE
*/
static cube curl2(mat psi) {
    cube u = grad(psi);
    cube temp = u;
    u.slice(0) = temp.slice(1);
    u.slice(1) = -temp.slice(0);
    return u;
}

/*
Created:        19 november 2013    by Jeroen Barnhoorn
Last changed:   10 december 2013    by Jeroen Barnhoorn

Replace dphi at the internal boundary with a 1-sided difference, on the
outside of the body.

STATUS: DONE
*/
static cube grad_bdry(mat phi, cube dphi, GridBdry gridbdry) {

    mat Gamma = gridbdry.Gamma;
    mat Gnormals = trans(gridbdry.n_grid);

    vec Gplus = zeros<vec>(Gamma.n_rows);
    vec Gminus = zeros<vec>(Gamma.n_rows);
    for(unsigned int i = 0; i < Gplus.n_rows; i++) {
        if(Gnormals(0, i) > 0) {
            Gplus(i) = Gamma(i, 0) + 1;
        } else {
            Gplus(i) = Gamma(i, 0);
        }
        if(Gnormals(0, i) > 0) {
            Gminus(i) = Gamma(i, 0);
        } else {
            Gminus(i) = Gamma(i, 0) - 1;
        }
    }

    for(unsigned int i = 0; i < Gplus.n_rows; i++) {
        dphi(Gamma(i, 0), Gamma(i, 1), 0) = phi(Gplus(i), Gamma(i, 1)) - phi(Gminus(i), Gamma(i, 1));
    }

    for(unsigned int i = 0; i < Gplus.n_rows; i++) {
        if(Gnormals(1, i) > 0) {
            Gplus(i) = Gamma(i, 1) + 1;
        } else {
            Gplus(i) = Gamma(i, 1);
        }
        if(Gnormals(1, i) > 0) {
            Gminus(i) = Gamma(i, 1);
        } else {
            Gminus(i) = Gamma(i, 1) - 1;
        }
    }

    for(unsigned int i = 0; i < Gplus.n_rows; i++) {
        dphi(Gamma(i, 0), Gamma(i, 1), 1) = phi(Gamma(i, 0), Gplus(i)) - phi(Gamma(i, 0), Gminus(i));
    }

    return dphi;
}

#endif // DIFFOPS_INCLUDED
