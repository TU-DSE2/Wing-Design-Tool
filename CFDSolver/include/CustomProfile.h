#ifndef CUSTOMPROFILE_H
#define CUSTOMPROFILE_H

#include "Profile.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <armadillo>
#include <math.h>
using namespace std;
using namespace arma;

/*
Created:        16 december 2013    by Jeroen Barnhoorn
Last updated:   16 december 2013    by Jeroen Barnhoorn

STATUS: WIP
*/

class CustomProfile:public Profile {
    public:
        CustomProfile();
        CustomProfile(mat, mat, mat, mat, mat, vec);
        void setSize(vec, vec, float, float, float, float);
    private:
        mat bdryPoints;pair<mat, int> open_file(string) ;
        mat scalecoodsdown(mat);
        void rotatecoods(float, mat&, int);
        void scalecoods(float, int, mat&);
        mat createoriginalnormals (int, mat, int);
        int pnpoly(int, mat, mat, float, float);
        mat ExactBoundary(mat, int);
        cube CreateGridBndryCube(mat, mat, int, int, int);
        mat CreateGridBndry(cube, int, int, int);
        mat FindOmega(mat, int , int, int);
};

#endif // CUSTOMPROFILE_H
