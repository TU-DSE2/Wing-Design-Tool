#ifndef VISCOUS_WALL_VORTICITY_FLUX_INCLUDED
#define VISCOUS_WALL_VORTICITY_FLUX_INCLUDED

#include <armadillo>
#include "math.h"

#include "ExactBdry.h"

using namespace arma;
using namespace std;

/*
Created:         4 december 2013    by Steve Brust
Last updated:   11 december 2013    by Jeroen Barnhoorn

STATUS: WIP
*/

static mat f_phi(cube inbound_X_cube, double d, float inbound_b){
    mat x = inbound_X_cube.slice(0);
    mat y = inbound_X_cube.slice(1);
    mat yswitch = zeros<mat>(y.n_rows, y.n_cols);
    mat nphi = zeros(x.n_rows, x.n_cols);
    //switch function, filling of phi matrix, and filling of deltaw matrix
    for (unsigned int i = 0; i < x.n_rows; i++){//y.n_elem fails b/c y is float
        for (unsigned int j = 0; j < x.n_cols; j++){
			if (y(i,j) >= 0){
                yswitch(i,j) = 1.0;
            }
			else {
                yswitch(i,j) = 0.0;
            }
            nphi(i, j) = 1/inbound_b * exp(-pow(y(i, j), 2)/pow(inbound_b, 2)) * (erf((d-x(i, j))/inbound_b) + erf(d+x(i,j))/inbound_b);
		}
    }
    return nphi % yswitch;
}

static mat viscous_wall_vorticity_flux(cube u, mat w, float deltat, float nu, ExactBdry exactbdry){
    //Add vorticity flux from no-slip wall to field.
    mat w_out = zeros<mat>(u.n_rows, u.n_cols);
    int M = exactbdry.npoints();
    float nup = nu;
    float deltatp = deltat*100;
    float b = sqrt(4*nup*deltatp);
    b = 1;  // Why?
    mat wloc = exactbdry.interp_tangent(u);
    int W = 2;  // Width of window
    cube X = zeros<cube>(2*W, 2*W, 2);
    mat deltaw = zeros<mat>(u.n_rows, u.n_cols);
    for(int m = 0; m < M; m++) {
        int i = exactbdry.locx(m);
        int j = exactbdry.locy(m);
        for(int n = 0; n < 2*W; n++) {
            X.slice(0).row(n) = trans(linspace(i-W+1, i+W, 2*W));
            X.slice(1).col(n) = linspace(j-W+1, j+W, 2*W);
        }
        mat phiX = f_phi(exactbdry.local_coords(X, m), exactbdry.ds(m), b);
        deltaw.submat(span(i-W+1, i+W), span(j-W+1, j+W)) += wloc(m)*phiX;
    }

    w_out = w + deltaw * 0.02;

    return w_out;
}

#endif //VISCOUS_WALL_VORTICITY_FLUX_INCLUDED
