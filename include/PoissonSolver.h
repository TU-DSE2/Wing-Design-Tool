#ifndef POISSON_INCLUDED
#define POISSON_INCLUDED

/*
Written by Jeroen Barnhoorn,
25 november 2013

Generates a Mat<int> Poisson matrix.

STATUS: DONE
*/
mat poissonMat(int size) {
    mat A(size, size, fill::zeros);
    int n = sqrt(size);
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            if(i == j) {
                A(i, j) = 4;
            }
            if(abs(i-j) == n) {
                A(i, j) = -1;
            }
            if(abs(i-j) == 1 && !(i%n == 0 && j == (i-1)) && !(j%n == 0 && i == (j-1))) {
                A(i, j) = -1;
            }
        }
    }
    return A;
}

/*
Written by Jeroen Barnhoorn,
27 november 2013

STATUS: DONE
*/
mat pSolve(mat w) {
    mat temppsi = zeros<mat>(w.n_rows, w.n_cols);
    mat tempw = trans(w);
    vec wvec = vectorise(tempw.submat(span(1, tempw.n_rows-2), span(1, tempw.n_cols-2)), 0);
    vec psivec = zeros<vec>(wvec.n_rows);
    mat pMat = poissonMat(wvec.n_rows);

    psivec = solve(pMat, wvec);

    mat psimat = reshape(psivec, w.n_rows-2, w.n_rows-2);
    temppsi.submat(span(1, temppsi.n_rows-2), span(1, temppsi.n_cols-2)) = psimat;

    return trans(temppsi);
}

#endif // POISSON_INCLUDED
