#include "CFDSolver.h"

CFDSolver::CFDSolver() {
    Nx = solverparameters.grid_size[0];
    Ny = solverparameters.grid_size[1];
    x = linspace<vec>(0, Nx-1, Nx);
    y = linspace<vec>(0, Ny-1, Ny);
    nu = solverparameters.viscosity;

    u = zeros<cube>(Nx, Ny, 2);
    w = zeros<mat>(Nx, Ny);
    ddw = zeros<mat>(Nx, Ny);
    psi = zeros<mat>(Nx, Ny);
    dphi = zeros<cube>(Nx, Ny, 2);
    g_force = zeros<mat>(Nx, Ny);

    profile.setSize(x, y, Nx/4, Ny/2, Nx/16, 120);
    gridbdry = profile.getGridBdry();
    exactbdry = profile.getExactBdry();

    g = zeros<vec>(gridbdry.npoints());

    Lmat = L_matrix();

    total_iterations = 0;
    iter = 0;
}

CFDSolver::CFDSolver(int Nxin, int Nyin) {
    Nx = Nxin;
    Ny = Nyin;
    x = linspace<vec>(0, Nx-1, Nx);
    y = linspace<vec>(0, Ny-1, Ny);
    nu = solverparameters.viscosity;

    u = zeros<cube>(Nx, Ny, 2);
    w = zeros<mat>(Nx, Ny);
    ddw = zeros<mat>(Nx, Ny);
    psi = zeros<mat>(Nx, Ny);
    dphi = zeros<cube>(Nx, Ny, 2);
    g_force = zeros<mat>(Nx, Ny);


    profile.setSize(x, y, Nx/4, Ny/2, Nx/16, 120);
    gridbdry = profile.getGridBdry();
    exactbdry = profile.getExactBdry();

    g = zeros<vec>(gridbdry.npoints());

    Lmat = L_matrix();

    total_iterations = 0;
    iter = 0;
}

rowvec CFDSolver::L_matvec(rowvec garg) {
    for(unsigned int i = 0; i < garg.n_cols; i++) {
        g_force(gridbdry.Gamma(i, 0), gridbdry.Gamma(i, 1)) = garg(i);
    }
    phi = pSolve(g_force);
    dphi = grad_bdry(phi, dphi, gridbdry);
    return gridbdry.interp_normal(dphi);
}

mat CFDSolver::L_matrix() {
    int M = gridbdry.npoints();
    Lmat = zeros<mat>(M, M);
    garg = zeros<rowvec>(M);
    for(int i = 0; i < M; i++) {
        garg = zeros<rowvec>(M);
        garg(i) = 1;
        L_matvec(garg);
        Lmat.col(i) = trans(L_matvec(garg));
    }
    return Lmat;
}

void CFDSolver::run(int iterations) {
    total_iterations = iter + iterations;
    while (iter < total_iterations) {
        psi = pSolve(-w);
        u = curl2(psi);

        u.slice(0) = u.slice(0) + solverparameters.mean_flow[0]*ones<mat>(u.slice(0).n_rows, u.slice(0).n_cols);
        u.slice(1) = u.slice(1) + solverparameters.mean_flow[1]*ones<mat>(u.slice(0).n_rows, u.slice(1).n_cols);

        if(solverparameters.boundaries) {
            b = -gridbdry.interp_normal(u);
            g = solve(Lmat, trans(b));
        }

        for(int i = 0; i < gridbdry.npoints(); i++) {
            g_force(gridbdry.Gamma(i, 0), gridbdry.Gamma(i, 1)) = g(i);
        }

        phi = pSolve(g_force);
        dphi = grad(phi);
        dphi = grad_bdry(phi, dphi, gridbdry);
        u += dphi;

        if(solverparameters.boundaries) {
            for(unsigned int i = 0; i < gridbdry.OmegaI.n_rows; i++) {
                u(gridbdry.OmegaI(i, 1), gridbdry.OmegaI(i, 0), 0) = 0;
                u(gridbdry.OmegaI(i, 1), gridbdry.OmegaI(i, 0), 1) = 0;
            }
            /*for(unsigned int i = 0; i < gridbdry.Gamma.n_rows; i++) {
                u(gridbdry.Gamma(i, 0), gridbdry.Gamma(i, 1), 0) = 0;
                u(gridbdry.Gamma(i, 0), gridbdry.Gamma(i, 1), 1) = 0;
            }*/
        }

        abs_u = abs(u);
        deltat = min(0.7/abs_u.max(), 1/(4.1*nu));
        for(int i = 0; i < solverparameters.n_advection_steps; i++) {
            w = interpolate_M4(w, u*deltat);
        }

        if(solverparameters.viscous) {
            if(solverparameters.boundaries) {
                w = viscous_wall_vorticity_flux(u, w, deltat, nu, exactbdry);
            }
            ddw = laplace(w);
            w += deltat*nu*ddw;
        }

        if(solverparameters.boundaries) {
            for(unsigned int i = 0; i < gridbdry.Gamma.n_rows; i++) {
                w(gridbdry.Gamma(i, 0), gridbdry.Gamma(i, 1)) *= 0;
            }
            for(unsigned int i = 0; i < gridbdry.OmegaI.n_rows; i++) {
                w(gridbdry.OmegaI(i, 1), gridbdry.OmegaI(i, 0)) *= 0;
            }
        }

        t += deltat;
        iter++;
    }
}

/*
  Generate a vortex at point (xcenter, ycenter).
*/

void CFDSolver::addVortex(int xcenter, int ycenter, double radius, double strength) {
    mat X = zeros<mat>(Nx, Nx);
    mat Y = zeros<mat>(Ny, Ny);
    for(unsigned int i = 0; i < X.n_rows; i++) {
        X.row(i) = trans(x);
    }
    for(unsigned int i = 0; i < Y.n_cols; i++) {
        Y.col(i) = y;
    }
    mat r = zeros<mat>(x.n_rows, y.n_rows);
    r = sqrt(square(trans(X)-xcenter) + square(trans(Y)-ycenter));
    mat temp = zeros<mat>(x.n_rows, y.n_rows);
    temp = 0.5*(cos(M_PI * r / radius) + 1)*strength;
    for(unsigned int i = 0; i < temp.n_cols; i++) {
        for(unsigned int j = 0; j < temp.n_rows; j++) {
            if(r(i, j) > radius) {
                temp(i, j) = 0;
            }
        }
    }
    w = w + temp;
}
