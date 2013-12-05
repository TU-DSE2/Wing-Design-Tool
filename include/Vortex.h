#ifndef VORTEX_INCLUDED
#define VORTEX_INCLUDED

using namespace std;
using namespace arma;

/*
Writtem by Jeroen Barnhoorn,
19 november 2013

Add a vortex to the vorticity field w.

STATUS: DONE
*/
mat add_vortex(mat w, vec x, vec y, int xcenter, int ycenter, double radius, double strength) {
    mat v = w;
    mat X = zeros<mat>(x.n_rows, x.n_rows);
    mat Y = zeros<mat>(y.n_rows, y.n_rows);
    for(unsigned int i = 0; i < X.n_rows; i++) {
        X.row(i) = trans(x);
    }
    for(unsigned int i = 0; i < Y.n_cols; i++) {
        Y.col(i) = y;
    }
    mat r = zeros<mat>(x.n_rows, y.n_rows);
    r = sqrt(square(trans(X)-xcenter) + square(trans(Y)-ycenter));
    mat temp = zeros<mat>(x.n_rows, y.n_rows);
    temp = 0.5*(cos(M_PI * r / 1.125) + 1)*strength;
    for(unsigned int i = 0; i < temp.n_cols; i++) {
        for(unsigned int j = 0; j < temp.n_rows; j++) {
            if(r(i, j) > radius) {
                temp(i, j) = 0;
            }
        }
    }
    v = v + temp;
    return v;
}


#endif // VORTEX_INCLUDED
