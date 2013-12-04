#ifndef FUNC_H_INCLUDED
#define FUNC_H_INCLUDED

/*
Written by: Jeroen Barnhoorn,
Created: 4 december 2013

STATUS: DONE
*/
mat nonzero(mat in) {
    umat temp = find(in);
    mat out = zeros<mat>(temp.n_rows, 2);
    out.col(0) = floor(conv_to<mat>::from(temp)/in.n_rows);
    for(unsigned int i = 0; i < temp.n_rows; i++) {
        out(i, 1) = temp(i)%in.n_rows;
    }
    return out;
}

#endif // FUNC_H_INCLUDED
