#ifndef DIFFOPS_INCLUDED
#define DIFFOPS_INCLUDED

#include <armadillo>

#include "GridBdry.h"

using namespace arma;

/*
Written by Jeroen Barnhoorn,
4 december 2013

STATUS: TODO
*/
mat laplace(mat w) {
    mat ddw = zeros<mat>(w.n_rows, w.n_cols);

    return ddw;
}

/*
Written by Jeroen Barnhoorn,
18 november 2013

STATUS: DONE
*/
cube grad(mat phi) {
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
Written by Jeroen Barnhoorn,
18 november 2013

Curl of vector => scalar

STATUS: DONE
*/
mat curl(cube u) {
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
Written by Jeroen Barnhoorn,
19 november 2013

Curl of scalar => vector

STATUS: DONE
*/
cube curl2(mat psi) {
    cube u = grad(psi);
    cube temp = u;
    u.slice(0) = temp.slice(1);
    u.slice(1) = -temp.slice(0);
    return u;
}

/*
Created: 19 november 2013 by Jeroen Barnhoorn
Last changed: 5 december 2013 by Jeroen Barnhoorn

Replace dphi at the internal boundary with a 1-sided difference, on the
outside of the body.

STATUS: DONE
*/
cube grad_bdry(mat phi, cube dphi, GridBdry gridbdry) {

    mat Gamma = gridbdry.Gamma;
    mat Gnormals = trans(gridbdry.n);

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
