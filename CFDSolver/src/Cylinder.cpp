#include "../include/Cylinder.h"

Cylinder::Cylinder(vec x, vec y, float xcenter, float ycenter, float radius, float M) {
    X = zeros<mat>(x.n_rows, x.n_rows);
    Y = zeros<mat>(y.n_rows, y.n_rows);
    for(unsigned int i = 0; i < X.n_rows; i++) {
        X.row(i) = trans(x);
    }
    for(unsigned int i = 0; i < Y.n_cols; i++) {
        Y.col(i) = y;
    }
    r = zeros<mat>(x.n_rows, y.n_rows);
    for(unsigned int i = 0; i < x.n_rows; i++) {
        for(unsigned int j = 0; j < y.n_rows; j++) {
            r(i, j) = sqrt(pow(Y(i, j) - xcenter, 2) + pow(X(i, j) - ycenter, 2));
        }
    }
    Oint_C = zeros<mat>(r.n_rows, r.n_cols);
    Gamma_C = zeros<mat>(r.n_rows, r.n_cols);
    for(unsigned int i = 0; i < Oint_C.n_rows; i++) {
        for(unsigned int j = 0; j < Oint_C.n_cols; j++) {
            if(r(i, j) <= (radius - 0.55)) {
                Oint_C(i, j) = 1;
            }
            if(r(i, j) <= (radius + 0.55)) {
                Gamma_C(i, j) = 1 - Oint_C(i, j);
            }
        }
    }
    Oint_N = sum(sum(Oint_C));
    Gamma_N = sum(sum(Gamma_C));
    interior = zeros<mat>(Oint_N, 2);
    Gamma_I = zeros<mat>(Gamma_N, 2);
    int n = 0;
    for(unsigned int i = 0; i < Gamma_C.n_rows; i++) {
        for(unsigned int j = 0; j < Gamma_C.n_cols; j++) {
            if(Gamma_C(i, j) != 0) {
                Gamma_I(n, 0) = i;
                Gamma_I(n, 1) = j;
                n++;
            }
        }
    }
    int m = 0;
    for(unsigned int i = 0; i < Oint_C.n_rows; i++) {
        for(unsigned int j = 0; j < Oint_C.n_cols; j++) {
            if(Oint_C(i, j) != 0) {
                interior(m, 0) = i;
                interior(m, 1) = j;
                m++;
            }
        }
    }
    n_grid_temp = zeros<mat>(Gamma_N, 2);
    n_grid = zeros<mat>(Gamma_N, 2);
    n_grid_temp.col(0) = Gamma_I.col(0) - xcenter;
    n_grid_temp.col(1) = Gamma_I.col(1) - ycenter;
    for(unsigned int i = 0; i < n_grid_temp.n_rows; i++) {
        n_grid(i, 0) = n_grid_temp(i, 0) / sqrt(pow(n_grid_temp(i, 0), 2)+pow(n_grid_temp(i, 1), 2));
        n_grid(i, 1) = n_grid_temp(i, 1) / sqrt(pow(n_grid_temp(i, 0), 2)+pow(n_grid_temp(i, 1), 2));
    }

    theta = linspace<vec>(0, 2*M_PI*(1-1/M), M);
    loc = zeros<mat>(M, 2);
    loc.col(0) = radius*cos(theta);
    loc.col(1) = radius*sin(theta);

    ex_ds = ones<mat>(M)*2*M_PI*radius/M;

    ex_n = zeros<mat>(loc.n_rows, loc.n_cols);
    ex_n.col(0) = loc.col(0)/radius;
    ex_n.col(1) = loc.col(1)/radius;

    ex_loc = zeros<mat>(loc.n_rows, loc.n_cols);
    ex_loc.col(0) = loc.col(0) + xcenter;
    ex_loc.col(1) = loc.col(1) + ycenter;
}

void Cylinder::setSize(vec x, vec y, float xcenter, float ycenter, float radius, float M) {
    X = zeros<mat>(x.n_rows, x.n_rows);
    Y = zeros<mat>(y.n_rows, y.n_rows);
    for(unsigned int i = 0; i < X.n_rows; i++) {
        X.row(i) = trans(x);
    }
    for(unsigned int i = 0; i < Y.n_cols; i++) {
        Y.col(i) = y;
    }
    r = zeros<mat>(x.n_rows, y.n_rows);
    for(unsigned int i = 0; i < x.n_rows; i++) {
        for(unsigned int j = 0; j < y.n_rows; j++) {
            r(i, j) = sqrt(pow(Y(i, j) - xcenter, 2) + pow(X(i, j) - ycenter, 2));
        }
    }
    Oint_C = zeros<mat>(r.n_rows, r.n_cols);
    Gamma_C = zeros<mat>(r.n_rows, r.n_cols);
    for(unsigned int i = 0; i < Oint_C.n_rows; i++) {
        for(unsigned int j = 0; j < Oint_C.n_cols; j++) {
            if(r(i, j) <= (radius - 0.55)) {
                Oint_C(i, j) = 1;
            }
            if(r(i, j) <= (radius + 0.55)) {
                Gamma_C(i, j) = 1 - Oint_C(i, j);
            }
        }
    }
    Oint_N = sum(sum(Oint_C));
    Gamma_N = sum(sum(Gamma_C));
    interior = zeros<mat>(Oint_N, 2);
    Gamma_I = zeros<mat>(Gamma_N, 2);
    int n = 0;
    for(unsigned int i = 0; i < Gamma_C.n_rows; i++) {
        for(unsigned int j = 0; j < Gamma_C.n_cols; j++) {
            if(Gamma_C(i, j) != 0) {
                Gamma_I(n, 0) = i;
                Gamma_I(n, 1) = j;
                n++;
            }
        }
    }
    int m = 0;
    for(unsigned int i = 0; i < Oint_C.n_rows; i++) {
        for(unsigned int j = 0; j < Oint_C.n_cols; j++) {
            if(Oint_C(i, j) != 0) {
                interior(m, 1) = i;
                interior(m, 0) = j;
                m++;
            }
        }
    }
    n_grid_temp = zeros<mat>(Gamma_N, 2);
    n_grid = zeros<mat>(Gamma_N, 2);
    n_grid_temp.col(0) = Gamma_I.col(0) - xcenter;
    n_grid_temp.col(1) = Gamma_I.col(1) - ycenter;
    for(unsigned int i = 0; i < n_grid_temp.n_rows; i++) {
        n_grid(i, 0) = n_grid_temp(i, 0) / sqrt(pow(n_grid_temp(i, 0), 2)+pow(n_grid_temp(i, 1), 2));
        n_grid(i, 1) = n_grid_temp(i, 1) / sqrt(pow(n_grid_temp(i, 0), 2)+pow(n_grid_temp(i, 1), 2));
    }

    theta = linspace<vec>(0, 2*M_PI*(1-1/M), M);
    loc = zeros<mat>(M, 2);
    loc.col(0) = radius*cos(theta);
    loc.col(1) = radius*sin(theta);

    ex_ds = ones<mat>(M)*2*M_PI*radius/M;

    ex_n = zeros<mat>(loc.n_rows, loc.n_cols);
    ex_n.col(0) = loc.col(0)/radius;
    ex_n.col(1) = loc.col(1)/radius;

    ex_loc = zeros<mat>(loc.n_rows, loc.n_cols);
    ex_loc.col(0) = loc.col(0) + xcenter;
    ex_loc.col(1) = loc.col(1) + ycenter;
}
