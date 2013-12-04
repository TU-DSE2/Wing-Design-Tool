#ifndef INTERPOLATE_INCLUDED
#define INTERPOLATE_INCLUDED

using namespace arma;

/*
Written by Jeroen Barnhoorn,
20 november 2013

Evaluate vorticity at grid-nodes given vorticity strength and coordinate
offset from nodes. Args: w - vorticity (Nx,Ny), xx node-offset (2,Nx,Ny).
Solution is returned in w.

Assumes that points are displaced by less than the mesh spacing in each
direction - i.e. |xx| < 1 for all elements of xx.  It may be interesting to
relax this condition to allow for bigger steps.  E.g. perhaps starting with
xa = np.fabs(xx[0])
xi = np.floor(xf)
xa -= xi
before computing the basis functions.
TODO: CPU intensive part of code - this will be the main overhead, together
with the Poisson solver.

STATUS: NEEDS DEBUGGING
*/
mat interpolate_M4(mat w, cube xx) {
    mat xa = abs(xx.slice(0));
    cube dx = zeros<cube>(w.n_rows, w.n_cols, 5);
    dx.slice(1) = (square(2 - (1 + xa))) % (1 - (1 + xa)) / 2;    // NOTE: % is element-wise multiplication
    dx.slice(2) = (3*pow(xa, 3) - 5*square(xa) + 2) / 2;
    dx.slice(3) = (3*pow((1 - xa), 3) - 5*square(1 - xa) + 2) / 2;
    dx.slice(4) = square(2 - (2 - xa)) % (1 - (2 - xa)) / 2;
    cube temp = dx;
    for(unsigned int i = 0; i < dx.n_rows; i++ ) {
        for(unsigned int j = 0; j < dx.n_cols; j++) {
            if(xx(i, j, 0) <= 0) {
                dx(i, j, 0) = temp(i, j, 4);
                dx(i, j, 1) = temp(i, j, 3);
                dx(i, j, 3) = temp(i, j, 1);
                dx(i, j, 4) = temp(i, j, 0);
            }
        }
    }

    mat ya = abs(xx.slice(1));
    cube dy = zeros<cube>(w.n_rows, w.n_cols, 5);
    dy.slice(1) = (square(2 - (1 + ya))) % (1 - (1 + ya)) / 2;
    dy.slice(2) = (3*pow(ya, 3) - 5*square(ya) + 2) / 2;
    dy.slice(3) = (3*pow((1 - ya), 3) - 5*square(1 - ya) + 2) / 2;
    dy.slice(4) = square(2 - (2 - ya)) % (1 - (2 - ya)) / 2;
    cube tempy = dy;

    for(unsigned int i = 0; i < dx.n_rows; i++ ) {
        for(unsigned int j = 0; j < dx.n_cols; j++) {
            if(xx(i, j, 1) <= 0) {
                dy(i, j, 0) = tempy(i, j, 4);
                dy(i, j, 1) = tempy(i, j, 3);
                dy(i, j, 3) = tempy(i, j, 1);
                dy(i, j, 4) = tempy(i, j, 0);
            }
        }
    }

    mat wnew = zeros<mat>(w.n_rows+4, w.n_cols+4);

    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            wnew.submat(span(i, i+w.n_rows-1), span(j, j+w.n_rows-1)) = wnew.submat(span(i, i+w.n_rows-1), span(j, j+w.n_rows-1)) + w%dx.slice(i)%dy.slice(j);
        }
    }

    return wnew.submat(span(2, wnew.n_rows-3), span(2, wnew.n_cols-3));
}

#endif // INTERPOLATE_INCLUDED
