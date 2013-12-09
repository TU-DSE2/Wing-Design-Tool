#ifndef VISCOUS_WALL_VORTICITY_FLUX_INCLUDED
#define VISCOUS_WALL_VORTICITY_FLUX_INCLUDED

#include <armadillo>
#include "math.h"

#include "ExactBdry.h"

using namespace arma;

/*
Written by Steve Brust
4 december 2013

STATUS: TODO
*/

/*
//phi is passed two vectors
//This function will go into the loop below. C++ does not allow functions within functions
float phi(vec x, float d){
	//vec y = zeros<vec>(x.n_elem);
	float(y) = x(1);//Takes a vector and creates two points, I think?
	x = x(0);
	//switch = np.where(y >= 0, 1., 0.)
	//Run through array y, if logic statement is true, set element to 1.0, if false set to 0.0

	//vec switch;
	for (int i = 0; i < 1; i++){//y.n_elem fails b/c y is float

	}
}*/

mat phi(cube x, mat d, float inbound_b){
	mat x = exactbdry.X.slice(0);
	mat y = exactbdry.X.slice(1);
	mat phi = zeros(x.n_rows, x.n_cols);
	//switch function, filling of phi matrix, and filling of deltaw matrix
	for (int i = 0; i < y.n_rows; i++){//y.n_elem fails b/c y is float
		for (int j = 0; j < y.n_cols; j++){
			if (y(i,j) >= 0){
				y(i,j) = 1.0;
			}
			else {
				y(i,j) = 0.0;
			}

			phi(i,j) = 1.0/inbound_b * exp(pow(-y(i,j), 2) / pow(inbound_b, 2))
				* (erfc((d(i) - x(i,j)) / inbound_b) + erfc((d(i) + x(i,j)) / inbound_b));

		}
	}

	return phi * y;
	//phiX = 1.0/b * exp(pow(-y,2)/pow(b,2)) * (erfc((d-x)/b) + erfc((d+x)/b))

}

void viscous_wall_vorticity_flux(cube u, mat w, float deltat, float nu, ExactBdry exactbdry){
	int M = exactbdry.npoints(); //number of point in exactbdry input
	float nup = nu;
	float deltatp = deltat*100;
	float b = sqrt(4.0*nup*deltatp);
	b = 1.0;//Why compute b but then set it to 1.0 immediately after?
	vec utau = exactbdry.interp_tangent(u);
	vec wloc = utau;
	int W = 2;
	cube X = zeros(2*W, 2*W, 2);
	mat deltaw = zeros(w.n_rows, w.n_cols);
	vec i;
	vec j;
	vec temp1; vec temp2;

	for (int m = 0; m < M; m++){
		i = exactbdry.locx.row(m);
		j = exactbdry.locy.row(m);

		temp1 = linspace<vec>(i - W + 1, i + W, 2*W);
		temp2 = linspace<vec>(j - W + 1, j + W, 2*W);

		//X[0], X[1] = np.meshgrid(np.linspace(i-W+1,i+W,2*W),np.linspace(j-W+1,j+W,2*W))
		//making meshgrid
		for (int i = 0; i < 2*W; i++){
			X.slice(0).row(i) = temp1;
			X.slice(1).col(i) = temp2;
		}

		mat phiX = phi(exactbdry.local_coords(X, m), exactbdry.ds.row(m), b);
		deltaw.submat(i-W+1, j-W+1, i+W+1, j+W+1) =
			deltaw.submat(i-W+1, j-W+1, i+W+1, j+W+1) + wloc.row(m) * phiX;

		//phi function below
		//Going to input exactbdry.local_coords and exactbdry.ds
		//phiX = phi(exactbdry.local_coords(X, m), exactbdry.ds[m];

		//the function, phi
		//vec y = zeros<vec>(x.n_elem);
			//float(y) = x(1);//Takes a vector and creates two points, I think?
			//x = x(0);

		//mat x = exactbdry.X.slice(0);
		//mat y = exactbdry.X.slice(1);

			//switch = np.where(y >= 0, 1., 0.)
			//Run through array y, if logic statement is true, set element to 1.0, if false set to 0.0

			//vec switch;
		/*
		for (int i = 0; i < y.n_rows; i++){//y.n_elem fails b/c y is float
			for (int j = 0; j < y.n_cols; j++){
				if (y(i,j) >= 0){
					y(i,j) = 1.0;
				}
				else {
					y(i,j) = 0.0;
				}
			}
		}

		for (int i = 0; i < ; i++){
			for (int j = 0; j < ; j++){
				mat phiX = 1.0/b * exp(pow(-y,2)/pow(b,2)) * (erfc((d-x)/b) + erfc((d+x)/b));

				deltaw.submat(i-W+1, j-W+1, i+W+1, j+W+1) =
						deltaw.submat(i-W+1, j-W+1, i+W+1, j+W+1) + wloc.row(m) * phiX

			}
		}


		//Now set value of phi
		//phi = 1.0/b * np.exp(-y**2/b**2) * (erf((d-x)/b) + erf((d+x)/b))
	*/
	}

	vec x1 = linspace(0, Nx - 1, Nx);
	vec y1 = linspace(0, Ny - 1, Ny);

	//MATfig statement omitted

	w = w + deltaw * 0.02;

	return w;

}

#endif //VISCOUS_WALL_VORTICITY_FLUX_INCLUDED
