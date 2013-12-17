#include <iostream>
#include <armadillo>
#include <paralution.hpp>
#include "CFDSolver/CFDSolver.h"
#include "CFDSolver/include/Cylinder.h"

using namespace std;
//using namespace arma;
//using namespace paralution;

int main() {
    CFDSolver cfdsolver(16, 16);
    cfdsolver.addVortex(cfdsolver.Nx*3/16, cfdsolver.Ny*1/2, (double)cfdsolver.Nx/8, 1);
    cfdsolver.run(1);
    //cout << cfdsolver.w;

    //arma::mat mat = arma::randu<arma::mat>(100,100);

    /*
    arma::vec val = arma::randu<arma::vec>(10);
    arma::ivec row = arma::linspace<arma::ivec>(0,9,10);//0-9
    arma::ivec col = arma::linspace<arma::ivec>(0,9,10);//0-9

    arma::vec * p_val = & val;
    arma::ivec * p_row = & row;
    arma::ivec * p_col = & col;


    cout << val << endl;
    cout << row << endl;
    cout << col << endl;
	*/
    /*
    cout << p_val << endl;//*p_val displays actual values
    cout << p_row << endl;
    cout << p_col << endl;
	*/


    arma::mat A = arma::zeros<arma::mat>(10,10);
        for (int i = 0; i < A.n_rows; i++){//makes identity matrix
        	for (int j = 0; j < A.n_cols; j++){
        		if (i == j){
        			A(i,j) = 1;
        		}
        	}
        }
    double* A_mem = A.memptr();

    arma::vec b = arma::ones<arma::vec>(10);
    double* b_mem = b.memptr();


    //////////////////////////////
    paralution::init_paralution();
    paralution::info_paralution();

    paralution::LocalMatrix<double> mat;
    mat.AllocateDENSE("testA", 10,10);
    /*
    arma::mat A = arma::zeros<arma::mat>(10,10);
    for (int i = 0; i < A.n_rows; i++){//makes identity matrix
    	for (int j = 0; j < A.n_cols; j++){
    		if (i == j){
    			A(i,j) = 1;
    		}
    	}
    }
    double* A_mem = A.memptr();
    */
    mat.SetDataPtrDENSE(&A_mem, "ones", 10, 10);
    mat.ConvertToCSR();
    mat.WriteFileMTX("test.mtx");

    paralution::LocalVector<double> vec;
    vec.SetValues(*b_mem);
    std::cout << "pos 1" << std::endl;
    vec.WriteFileASCII("test.txt");
    std::cout << "pos 2" << std::endl;

    //vec.AllocateDense("testb",10);





    //paralution::LocalMatrix<double> mat1;
    //mat1.paralution::LocalMatrix<double>::CopyFromCOO(*p_row,*p_col,*p_val);

    //paralution::LocalMatrix<double> para_mat;
    //mat.ReadFileMTX('gr_30_30.mtx');
    //para_mat.SetDataPtr(&arma_mat, "vector",200);
    /*
    paralution::LocalMatrix<double> mat;
    mat.AllocateCSR("My CSR Matrix", 5, 5, 5);

    mat.CopyFromCOO(row,col,val);
    //mat.Transpose();
	*/

    for (int i = 0 ; i <= 7; i++){
    	std::cout << i << endl;
    }


    paralution::stop_paralution();
    //////////////////////////////
    std::cout << "done" << std::endl;

    return 0;
}
