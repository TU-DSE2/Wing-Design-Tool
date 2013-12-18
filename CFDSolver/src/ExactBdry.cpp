/*
 * ExactBdry.cpp
 *
 *  Created on: Dec 5, 2013
 *      Author: Steve Brust
 *  Last updated: Dec 18, 2013 (S. Brust)
 */

#include "../include/ExactBdry.h"

ExactBdry::ExactBdry(mat inbound_loc, mat inbound_ds, mat inbound_n) {
	loc = inbound_loc;
	ds = inbound_ds;
	n = inbound_n;

	tau = zeros(n.n_rows,2);
	tau.col(0) = n.col(1);
	tau.col(1) = -n.col(0);

	mat loci = floor(loc);

	locx = loci.col(0);
	locy = loci.col(1);

	loco = loc - loci;

	M = loc.n_rows;

    basis = zeros(4, M);
	for (int i = 0; i < M; i++){
		basis(0,i) = (1 - loco(i,0)) * (1 - loco(i,1));
		basis(1,i) = loco(i,0) * (1 - loco(i,1));
		basis(2,i) = loco(i,0) * loco(i,1);
		basis(3,i) = (1 - loco(i,0)) * loco(i,1);
	}
}

int ExactBdry::npoints(){
	return M;
}

mat ExactBdry::interp(cube inbound_u){
	u = inbound_u;
	mat temp = zeros(2,M);

    for(int i = 0; i < M; i++){
        for(unsigned int j = 0; j < u.n_slices; j++){

			temp(j,i) = temp(j,i)
					+ basis(0,i)*u(locx(i),locy(i),j)
					+ basis(1,i)*u(locx(i) + 1,locy(i),j)
					+ basis(2,i)*u(locx(i) + 1,locy(i) + 1,j)
					+ basis(3,i)*u(locx(i),locy(i) + 1,j);
		}
	}

	return temp;
}

vec ExactBdry::interp_tangent(cube inbound_u){
	u = inbound_u;
	mat iu = interp(u);

	vec temp_interp_tangent = zeros<vec>(tau.n_rows);  //Supposed to be equal to M
	if (iu.n_rows == tau.n_cols && tau.n_rows == iu.n_cols){
        for (unsigned int i = 0; i < tau.n_rows; i++){//will run to M
            for (unsigned int j = 0; j < tau.n_cols; j++){//will run to 2
				temp_interp_tangent(i) = temp_interp_tangent(i) + (iu(j,i) * tau(i,j));
			}
		}
	}

	else {
        //cout << "The dimensions of tau and u were not compatible" << endl;
	}

	return temp_interp_tangent;
}

vec ExactBdry::interp_normal(cube inbound_u){
	u = inbound_u;
	mat iu = interp(u);

	vec temp_interp_normal = zeros<vec>(n.n_rows);  //Supposed to be equal to M
	if (iu.n_rows == n.n_cols && n.n_rows == iu.n_cols){
        for (unsigned int i = 0; i < n.n_rows; i++){//will run to M
            for (unsigned int j = 0; j < n.n_cols; j++){//will run to 2
				temp_interp_normal(i) = temp_interp_normal(i) + (iu(j,i) * n(i,j));
			}
		}
	}

	else {
        //cout << "The dimensions of n and u were not compatible" << endl;
	}

	return temp_interp_normal;
}

cube ExactBdry::local_coords(cube inbound_x, int inbound_m){
	x = inbound_x;//public vars equal to the input values
    m = inbound_m;

	mat ntau = zeros(tau.n_cols,n.n_cols);//Usually a 2x2 matrix
	ntau.row(0) = tau.row(m);
	ntau.row(1) = n.row(m);

	//setting up mats for x_transpose of cube and temps
	cube x_trans = zeros(x.n_rows, x.n_slices, x.n_cols);//arma does not have cube.t() func
	cube temp1 = zeros(x.n_rows, x.n_slices, x.n_cols);//x.T-self.loc[m] term
	cube temp2 = zeros(temp1.n_rows, temp1.n_slices, temp1.n_cols);//(x.T-self.loc[m]).T term

	//used to evaluate x transpose (of cube) and the x.T-self.loc[m] term
    for (unsigned int i = 0; i < x.n_rows; i++){
        for (unsigned int j = 0; j < x.n_slices; j++ ){
            for (unsigned int k = 0; k < x.n_cols; k++){
				x_trans(i,j,k) = x(i,k,j);
				temp1(i,j,k) = x_trans(i,j,k) - loc(m,j);
				//temp1(i,j,k) = x_trans(i,j,k) - loc(0,j);
            }
		}
	}

	//filling the temp2 matrix and finding the transpose of x.T-self.loc[m]
    for (unsigned int i = 0; i < temp1.n_rows; i++){
        for (unsigned int j = 0; j < temp1.n_slices; j++ ){
            for (unsigned int k = 0; k < temp1.n_cols; k++){
				temp2(i,j,k) = temp1(i,k,j);
			}
		}
	}

	cube result = zeros(4,4,2);//initialize matrix to fill with result

	//the following loop is the .py line np.einsum("ijk,li->ljk",(x.T-self.loc[m].T ...
	//A bit of hardcoding in here ntau(0,0) and ntau(0,1)
    for(unsigned int i = 0; i < result.slice(0).n_rows; i++){
        for(unsigned int j = 0; j < result.slice(0).n_cols; j++){
			result(i,j,0) = temp2(i,j,0)*ntau(0,0) + temp2(i,j,1)*ntau(0,1);
			result(i,j,1) = temp2(i,j,0)*ntau(1,0) + temp2(i,j,1)*ntau(1,1);
		}//result(i,j,0) = temp2(i,j,0)*ntau(0,0) + temp2(i,j,1)*ntau(1,0);
	}
	return result;
}

ExactBdry::~ExactBdry() {
	// TODO Auto-generated destructor stub
}

