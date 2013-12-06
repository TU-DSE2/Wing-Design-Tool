<<<<<<< HEAD
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

	ExactBdry(mat inbound_loc, mat inbound_ds, mat inbound_n);

	int npoints();
	mat interp(cube inbound_u);
	vec interp_tangent(cube inbound_u);
	vec interp_normal(cube inound_u);
	cube local_coords(cube inbound_x, int inbound_m);

	virtual ~ExactBdry();
};
#endif /* EXACTBDRY_H_ */
=======
#ifndef EXACTBDRY_H
#define EXACTBDRY_H

#include <armadillo>

using namespace arma;


class ExactBdry {
    public:
        ExactBdry(mat, vec, mat);
};

#endif // EXACTBDRY_H
>>>>>>> 57cf9c674d580d494b6bf5aba197d5d778cdd933
