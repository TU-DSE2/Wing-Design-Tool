#ifndef VISCOUS_WALL_VORTICITY_FLUX_INCLUDED
#define VISCOUS_WALL_VORTICITY_FLUX_INCLUDED

#include <armadillo>
#include "math.h"

#include "ExactBdry.h"

using namespace arma;

/*
Written by Steve Brust
4 december 2013

STATUS: Done; in need of validation and implementation
*/

static mat f_phi(cube inbound_X_cube, double d, float inbound_b){
    mat x = inbound_X_cube.slice(0);
    mat y = inbound_X_cube.slice(1);
    mat nphi = zeros(x.n_rows, x.n_cols);
    //switch function, filling of phi matrix, and filling of deltaw matrix
    for (unsigned int i = 0; i < x.n_rows; i++){//y.n_elem fails b/c y is float
        for (unsigned int j = 0; j < x.n_cols; j++){
			if (y(i,j) >= 0){
				y(i,j) = 1.0;
            }
			else {
				y(i,j) = 0.0;
            }
            nphi(i,j) = 1.0/inbound_b * exp(pow(-y(i,j), 2) / pow(inbound_b, 2))
				* (erfc((d - x(i,j)) / inbound_b) + erfc((d + x(i,j)) / inbound_b));
		}
    }

    return nphi * y;
}

static mat viscous_wall_vorticity_flux(cube u, mat w, float deltat, float nu, ExactBdry exactbdry){
	int M = exactbdry.npoints(); //number of point in exactbdry input
    float nup = nu;
	float deltatp = deltat*100;
    float b = sqrt(4.0*nup*deltatp);
    b = 1.0;//Why compute b but then set it to 1.0 immediately after?
    vec utau = exactbdry.interp_tangent(u);
	vec wloc = utau;
	int W = 2;
	cube X = zeros(2*W, 2*W, 2);
    mat deltaw = zeros(w.n_rows, w.n_cols);//9x9
	int i;
    int j;
    vec temp1; vec temp2;
    for (int m = 0; m < M; m++){
		i = exactbdry.locx(m);
		j = exactbdry.locy(m);

		vec temp1 = linspace<vec>(i - W + 1, i + W, 2*W);
		vec temp2 = linspace<vec>(j - W + 1, j + W, 2*W);

		//X[0], X[1] = np.meshgrid(np.linspace(i-W+1,i+W,2*W),np.linspace(j-W+1,j+W,2*W))
        //making meshgrid
		for (int k = 0; k < 2*W; k++){
            X.slice(0).row(k) = trans(temp1);
            X.slice(1).col(k) = temp2;
        }

		//function call to above f_phi function
        mat phiX = f_phi(exactbdry.local_coords(X, m), exactbdry.ds(m), b);

        //The following two for loops is simply the line
        //deltaw[i-W+1:i+W+1,j-W+1:j+W+1] += wloc[m] * phiX
        mat temp_mat = zeros(phiX.n_rows, phiX.n_cols);
        for (int k = 0; k < (i+W+1) - (i-W+1); k++){//k should run from 0 to 3
        	for (int l = 0; l < (j+W+1) - (j-W+1); l++){//j should run from 0 to 3
        		temp_mat(k,l) = deltaw(i-W+1 + k, j-W+1 + l);
        	}
        }

        temp_mat = temp_mat + (wloc(m) * phiX);

        for (int k = 0; k < (i+W+1) - (i-W+1); k++){//k should run from 0 to 3
            for (int l = 0; l < (j+W+1) - (j-W+1); l++){//j should run from 0 to 3
            	deltaw(i-W+1 + k, j-W+1 + l) = temp_mat(k,l);
            }
        }
    }

    w = w + deltaw * 0.02;
    return w;

}

#endif //VISCOUS_WALL_VORTICITY_FLUX_INCLUDED
