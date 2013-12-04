#include "../include/Rod.h"
#include "../include/Func.h"
#include "../include/GridBdry.h"

/*
Written by: Jeroen Barnhoorn,
Created: 4 december 2013

STATUS: WIP
*/
Rod::Rod(vec x, vec y) {
    //Define dimensions
    int Nx = x.n_rows;
    int Ny = y.n_rows;
    int Nx8 = floor(Nx/8);
    int Nx2 = floor(Nx/2);
    int Ny2 = floor(Ny/2);
    int thick = 2;
    int length = Nx2 - Nx8 + 1;

    //Define boundary corners
    rowvec corner_topleft(2);
    rowvec corner_topright(2);
    rowvec corner_botleft(2);
    rowvec corner_botright(2);
    corner_topleft << Nx8 << Ny2+thick;
    corner_topright << Nx2 << Ny2+thick;
    corner_botright << Nx2 << Ny2-thick;
    corner_botleft << Nx8 << Ny2-thick;

    //Define boundary sides
    mat side_top(corner_topright(0)-corner_topleft(0)-1, 2);
    mat side_bot(corner_botright(0)-corner_botleft(0)-1, 2);
    mat side_left(corner_topleft(1)-corner_botleft(1)-1, 2);
    mat side_right(corner_topright(1)-corner_botright(1)-1, 2);
    side_top.col(0) = linspace<vec>(Nx8+1, Nx2-1, length-2);
    side_top.col(1) = ones<vec>(length-2)*Ny2+thick;
    side_bot.col(0) = linspace<vec>(Nx2-1, Nx8+1, length-2);
    side_bot.col(1) = ones<vec>(length-2)*Ny2-thick;
    side_left.col(0) = ones<vec>(2*thick-1)*Nx8;
    side_left.col(1) = linspace<vec>(Ny2-thick+1, Ny2+thick-1, thick+1);
    side_right.col(0) = ones<vec>(2*thick-1)*Nx2;
    side_right.col(1) = linspace<vec>(Ny2+thick-1, Ny2-thick+1, thick+1);

    //Define normal vectors
    mat normal_top = zeros<mat>(length-2, 2);
    normal_top.col(1) = ones<vec>(length-2);
    mat normal_bot = zeros<mat>(length-2, 2);
    normal_bot.col(1) = -ones<vec>(length-2);
    mat normal_left = zeros<mat>(2*thick-1, 2);
    normal_left.col(0) = -ones<vec>(2*thick-1);
    mat normal_right = zeros<mat>(2*thick-1, 2);
    normal_right.col(0) = ones<vec>(2*thick-1);
    rowvec normal_topleft(2);
    rowvec normal_topright(2);
    rowvec normal_botleft(2);
    rowvec normal_botright(2);
    normal_topleft << -1/sqrt(2) << 1/sqrt(2);
    normal_topright << 1/sqrt(2) << 1/sqrt(2);
    normal_botleft << 1/sqrt(2) << -1/sqrt(2);
    normal_botright << -1/sqrt(2) << -1/sqrt(2);

    //Combine boundary points
    Gamma_I = zeros<mat>(side_top.n_rows+side_bot.n_rows+side_left.n_rows+side_right.n_rows+4, 2);
    Gamma_I.row(0) = corner_topleft;
    Gamma_I.submat(span(1, side_top.n_rows), span(0, 1)) = side_top;
    Gamma_I.row(1+side_top.n_rows) = corner_topright;
    Gamma_I.submat(span(2+side_top.n_rows, 1+side_top.n_rows+side_right.n_rows), span(0, 1)) = side_right;
    Gamma_I.row(2+side_top.n_rows+side_right.n_rows) = corner_botright;
    Gamma_I.submat(span(3+side_top.n_rows+side_right.n_rows, 2+side_top.n_rows+side_right.n_rows+side_bot.n_rows), span(0, 1)) = side_bot;
    Gamma_I.row(3+side_top.n_rows+side_right.n_rows+side_bot.n_rows) = corner_botleft;
    Gamma_I.submat(span(4+side_top.n_rows+side_right.n_rows+side_bot.n_rows, 3+side_top.n_rows+side_right.n_rows+side_bot.n_rows+side_left.n_rows), span(0, 1)) = side_left;

    //Combine normal vectors
    n = zeros<mat>(normal_top.n_rows+normal_bot.n_rows+normal_left.n_rows+normal_right.n_rows+4, 2);
    n.row(0) = normal_topleft;
    n.submat(span(1, normal_top.n_rows), span(0, 1)) = normal_top;
    n.row(1+normal_top.n_rows) = normal_topright;
    n.submat(span(2+normal_top.n_rows, 1+normal_top.n_rows+normal_right.n_rows), span(0, 1)) = normal_right;
    n.row(2+normal_top.n_rows+normal_right.n_rows) = normal_botright;
    n.submat(span(3+normal_top.n_rows+normal_right.n_rows, 2+normal_top.n_rows+normal_right.n_rows+normal_bot.n_rows), span(0, 1)) = normal_bot;
    n.row(3+normal_top.n_rows+normal_right.n_rows+normal_bot.n_rows) = normal_botleft;
    n.submat(span(4+normal_top.n_rows+normal_right.n_rows+normal_bot.n_rows, 3+normal_top.n_rows+normal_right.n_rows+normal_bot.n_rows+normal_left.n_rows), span(0, 1)) = normal_left;

    mat interiortemp = zeros<mat>(Nx, Ny);
    interiortemp.submat(span(Nx8+1, Nx2-1), span(Ny2-thick+1, Ny2+thick-1)) = ones<mat>(Nx2-Nx8-1, 2*thick-1);

    interior = nonzero(interiortemp);
    /*
                                        # Exact bdry points
    exact_top   = np.vstack((np.arange(Nx8,Nx2)+   .5, np.ones(length-1)*Ny2+thick)).T
    exact_bot   = np.vstack((np.arange(Nx2,Nx8,-1)-.5, np.ones(length-1)*Ny2-thick)).T
    exact_left  = np.vstack((np.ones(2*thick)*Nx8,    np.arange(Ny2-thick,Ny2+thick)   +.5)).T
    exact_right = np.vstack((np.ones(2*thick)*Nx2,    np.arange(Ny2+thick,Ny2-thick,-1)-.5)).T
    exnormal_top = np.array([0., 1.]) * np.ones((length-1,1))
    exnormal_bot   = np.array([0.,-1.]) * np.ones((length-1,1))
    exnormal_left  = np.array([-1.,0.]) * np.ones((2*thick,1))
    exnormal_right = np.array([ 1.,0.]) * np.ones((2*thick,1))

    ex_loc = np.vstack((exact_top, exact_right, exact_bot, exact_left))
    ex_n = np.vstack((exnormal_top, exnormal_right, exnormal_bot, exnormal_left))
    ex_ds  = np.ones(ex_n.shape[0]) * 1.

    return (solver.GridBdry(Gamma_I, np.nonzero(interior), n),
            solver.ExactBdry(ex_loc, ex_ds, ex_n))*/
}

GridBdry Rod::getGridBdry() {
    GridBdry gridbdry(Gamma_I, interior, n);
    return gridbdry;
}
