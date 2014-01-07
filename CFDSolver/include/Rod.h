#ifndef ROD_H
#define ROD_H

#include <armadillo>
#include "GridBdry.h"
#include "ExactBdry.h"

using namespace arma;

/*
Created:        4 december  by Jeroen Barnhoorn
Last updated:   10 december by Jeroen Barnhoorn

STATUS: WIP
*/
class Rod {
    public:
        vec x;
        vec y;
        float float1;
        float float2;
        float float3;
        float float4;

        Rod(vec, vec);
        Rod() { }

        void setSize(vec, vec, float, float, float, float);
        GridBdry getGridBdry();
        ExactBdry getExactBdry();
    protected:
    private:
        mat Gamma_I;
        mat n;
        mat interior;
        mat ex_loc;
        mat ex_n;
        vec ex_ds;
};

#endif // ROD_H
