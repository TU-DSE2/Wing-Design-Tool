#ifndef CYLINDER_H
#define CYLINDER_H

#include "Profile.h"

/*
Created:        13 december 2013    by Jeroen Barnhoorn
Last updated:   13 december 2013    by Jeroen Barnhoorn

STATUS: DONE
*/

class Cylinder:public Profile {
    public:
        Cylinder() { }
        Cylinder(vec, vec, float, float, float, float);
        void setSize(vec, vec, float, float, float, float);
    private:
        vec x;
        vec y;
        mat X;
        mat Y;
        mat r;
        mat Oint_C;
        mat Gamma_C;
        int Oint_N;
        int Gamma_N;
        mat n_grid_temp;
        vec theta;
        mat loc;
};

#endif // CYLINDER_H
