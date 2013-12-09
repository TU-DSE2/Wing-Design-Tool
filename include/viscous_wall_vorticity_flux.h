#ifndef VISCOUS_WALL_VORTICITY_FLUX_INCLUDED
#define VISCOUS_WALL_VORTICITY_FLUX_INCLUDED

#include <armadillo>

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

		temp1 = linspace(i - W + 1, i + W, 2*W);
		temp2 = linspace(j - W + 1, j + W, 2*W);

		//X[0], X[1] = np.meshgrid(np.linspace(i-W+1,i+W,2*W),np.linspace(j-W+1,j+W,2*W))

		//meshgrid still needs to happen

		//X.slice(0) = ;
		//X.slice(1) = ;
	}


}

#endif //VISCOUS_WALL_VORTICITY_FLUX_INCLUDED
