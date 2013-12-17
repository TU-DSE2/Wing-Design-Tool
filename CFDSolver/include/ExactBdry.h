/*
 * ExactBdry.h
 *
 *  Created on: Dec 5, 2013
 *      Author: Steve Brust
 */

#ifndef EXACTBDRY_H_
#define EXACTBDRY_H_

#include <armadillo>
using namespace arma;

class ExactBdry {
public:

	mat loc;
	mat ds;
	mat n;

	mat tau;
	vec locx;
	vec locy;
	mat loco;
	int M;
	mat basis;

	cube u;

	cube x;
	int m;

    ExactBdry() { }
    ExactBdry(mat, mat, mat);

	int npoints();
	mat interp(cube inbound_u);
	vec interp_tangent(cube inbound_u);
	vec interp_normal(cube inound_u);
	cube local_coords(cube inbound_x, int inbound_m);

	virtual ~ExactBdry();
};
#endif /* EXACTBDRY_H_ */
