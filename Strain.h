#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <complex>
#include <math.h>
#include <random>

#define MKL_Complex16 std::complex<double>
#define MKL_Complex8  std::complex<float>
#include <mkl.h>
#include "mkl_spblas.h"

#include "IO_data.h"


using namespace std;

#ifndef Strain_H

class Strain
{

public:
    Strain();
    ~Strain();

    int uniax_flag, shear_flag;
    double poisson, uniax_val, shear_val;

    void get_input_strain_data();
    
    void add_strain(int n,
                    int n_x,
                    complex<double> *h_2_u,
                    complex<double> *h_2_d,
                    complex<double> *h_3_u,
                    complex<double> *h_3_d,
                    complex<double> *x,
                    complex<double> *y);
                    
    void uniax_modify_hamiltonian(int n,
                                  complex<double> *h_2_u,
                                  complex<double> *h_2_d,
                                  complex<double> *h_3_u,
                                  complex<double> *h_3_d);
								  
    void uniax_modify_x_y(int n,
                          complex<double> *x,
                          complex<double> *y);
						  
    void shear_modify_hamiltonian(int n,
                                  int n_x,
                                  complex<double> *h_2_u,
                                  complex<double> *h_2_d,
                                  complex<double> *h_3_u,
                                  complex<double> *h_3_d,
                                  complex<double> *x,
                                  complex<double> *y);
								  
    void shear_modify_x_y(int n,
                          complex<double> *x,
                          complex<double> *y);
						  
    double get_distance_between_two_points(double x1,
                                           double y1,
                                           double x2,
                                           double y2); 

private:
    IO_data io_data;
};

#define Strain_H
#endif // !Strain_H




