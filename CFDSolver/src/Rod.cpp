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
    normal_botleft << -1/sqrt(2) << -1/sqrt(2);
    normal_botright << 1/sqrt(2) << -1/sqrt(2);

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

    //Define exact boundary sides
    mat exact_top(length-1, 2);
    exact_top.col(1) = ones<mat>(length-1)*Ny2+thick;
    exact_top.col(0) = linspace<vec>(Nx8+0.5, Nx2-0.5, length-1);
    mat exact_bot(length-1, 2);
    exact_bot.col(1) = ones<mat>(length-1)*Ny2-thick;
    exact_bot.col(0) = linspace<vec>(Nx2-0.5, Nx8+0.5, length-1);
    mat exact_left(2*thick, 2);
    exact_left.col(0) = ones<mat>(2*thick)*Nx8;
    exact_left.col(1) = linspace<vec>(Ny2-thick+0.5, Ny2+thick-0.5, 2*thick);
    mat exact_right(2*thick, 2);
    exact_right.col(0) = ones<mat>(2*thick)*Nx2;
    exact_right.col(1) = linspace<vec>(Ny2+thick-0.5, Ny2-thick+0.5, 2*thick);

    //Define exact normal vectors
    mat exnormal_top = zeros<mat>(length-1, 2);
    exnormal_top.col(1) = ones<vec>(length-1);
    mat exnormal_bot = zeros<mat>(length-1, 2);
    exnormal_bot.col(1) = -ones<vec>(length-1, 1);
    mat exnormal_left = zeros<mat>(2*thick, 2);
    exnormal_left.col(0) = -ones<vec>(2*thick);
    mat exnormal_right = zeros<mat>(2*thick, 2);
    exnormal_right.col(0) = ones<vec>(2*thick);

    //Combine exact boundary sides
    ex_loc = zeros<mat>(exact_top.n_rows+exact_bot.n_rows+exact_left.n_rows+exact_right.n_rows, 2);
    ex_loc.submat(span(0, exact_top.n_rows-1), span(0, 1)) = exact_top;
    ex_loc.submat(span(exact_top.n_rows, exact_top.n_rows+exact_right.n_rows-1), span(0, 1)) = exact_right;
    ex_loc.submat(span(exact_top.n_rows+exact_right.n_rows, exact_top.n_rows+exact_right.n_rows+exact_bot.n_rows-1), span(0, 1)) = exact_bot;
    ex_loc.submat(span(exact_top.n_rows+exact_right.n_rows+exact_bot.n_rows, exact_top.n_rows+exact_right.n_rows+exact_bot.n_rows+exact_left.n_rows-1), span(0, 1)) = exact_left;

    //Combine exact normal vectors
    ex_n = zeros<mat>(exnormal_top.n_rows+exnormal_bot.n_rows+exnormal_left.n_rows+exnormal_right.n_rows, 2);
    ex_n.submat(span(0, exnormal_top.n_rows-1), span(0, 1)) = exnormal_top;
    ex_n.submat(span(exnormal_top.n_rows, exnormal_top.n_rows+exnormal_right.n_rows-1), span(0, 1)) = exnormal_right;
    ex_n.submat(span(exnormal_top.n_rows+exnormal_right.n_rows, exnormal_top.n_rows+exnormal_right.n_rows+exnormal_bot.n_rows-1), span(0, 1)) = exnormal_bot;
    ex_n.submat(span(exnormal_top.n_rows+exnormal_right.n_rows+exnormal_bot.n_rows, exnormal_top.n_rows+exnormal_right.n_rows+exnormal_bot.n_rows+exnormal_left.n_rows-1), span(0, 1)) = exnormal_left;

    ex_ds = ones<vec>(ex_n.n_rows);
}

GridBdry Rod::getGridBdry() {
    GridBdry gridbdry(Gamma_I, interior, n);
    return gridbdry;
}

ExactBdry Rod::getExactBdry() {
    ExactBdry exactbdry(ex_loc, ex_ds, ex_n);
    return exactbdry;
}

void Rod::setSize(vec x_in, vec y_in, float float1, float float2, float float3, float float4) {
    x = x_in;
    y = y_in;

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
    normal_botleft << -1/sqrt(2) << -1/sqrt(2);
    normal_botright << 1/sqrt(2) << -1/sqrt(2);

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

    //Define exact boundary sides
    mat exact_top(length-1, 2);
    exact_top.col(1) = ones<mat>(length-1)*Ny2+thick;
    exact_top.col(0) = linspace<vec>(Nx8+0.5, Nx2-0.5, length-1);
    mat exact_bot(length-1, 2);
    exact_bot.col(1) = ones<mat>(length-1)*Ny2-thick;
    exact_bot.col(0) = linspace<vec>(Nx2-0.5, Nx8+0.5, length-1);
    mat exact_left(2*thick, 2);
    exact_left.col(0) = ones<mat>(2*thick)*Nx8;
    exact_left.col(1) = linspace<vec>(Ny2-thick+0.5, Ny2+thick-0.5, 2*thick);
    mat exact_right(2*thick, 2);
    exact_right.col(0) = ones<mat>(2*thick)*Nx2;
    exact_right.col(1) = linspace<vec>(Ny2+thick-0.5, Ny2-thick+0.5, 2*thick);

    //Define exact normal vectors
    mat exnormal_top = zeros<mat>(length-1, 2);
    exnormal_top.col(1) = ones<vec>(length-1);
    mat exnormal_bot = zeros<mat>(length-1, 2);
    exnormal_bot.col(1) = -ones<vec>(length-1, 1);
    mat exnormal_left = zeros<mat>(2*thick, 2);
    exnormal_left.col(0) = -ones<vec>(2*thick);
    mat exnormal_right = zeros<mat>(2*thick, 2);
    exnormal_right.col(0) = ones<vec>(2*thick);

    //Combine exact boundary sides
    ex_loc = zeros<mat>(exact_top.n_rows+exact_bot.n_rows+exact_left.n_rows+exact_right.n_rows, 2);
    ex_loc.submat(span(0, exact_top.n_rows-1), span(0, 1)) = exact_top;
    ex_loc.submat(span(exact_top.n_rows, exact_top.n_rows+exact_right.n_rows-1), span(0, 1)) = exact_right;
    ex_loc.submat(span(exact_top.n_rows+exact_right.n_rows, exact_top.n_rows+exact_right.n_rows+exact_bot.n_rows-1), span(0, 1)) = exact_bot;
    ex_loc.submat(span(exact_top.n_rows+exact_right.n_rows+exact_bot.n_rows, exact_top.n_rows+exact_right.n_rows+exact_bot.n_rows+exact_left.n_rows-1), span(0, 1)) = exact_left;

    //Combine exact normal vectors
    ex_n = zeros<mat>(exnormal_top.n_rows+exnormal_bot.n_rows+exnormal_left.n_rows+exnormal_right.n_rows, 2);
    ex_n.submat(span(0, exnormal_top.n_rows-1), span(0, 1)) = exnormal_top;
    ex_n.submat(span(exnormal_top.n_rows, exnormal_top.n_rows+exnormal_right.n_rows-1), span(0, 1)) = exnormal_right;
    ex_n.submat(span(exnormal_top.n_rows+exnormal_right.n_rows, exnormal_top.n_rows+exnormal_right.n_rows+exnormal_bot.n_rows-1), span(0, 1)) = exnormal_bot;
    ex_n.submat(span(exnormal_top.n_rows+exnormal_right.n_rows+exnormal_bot.n_rows, exnormal_top.n_rows+exnormal_right.n_rows+exnormal_bot.n_rows+exnormal_left.n_rows-1), span(0, 1)) = exnormal_left;

    ex_ds = ones<vec>(ex_n.n_rows);
}
