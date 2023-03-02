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
#include "Hamiltonian.h"


using namespace std;

#ifndef Psy_H

class Psy
{

public:
    Psy(Hamiltonian &h);
    ~Psy();

    complex<double> *psy_null;
    complex<double> *psy_null_dos;
    complex<double> *psy_x_null;

    complex<double> *psy;
    complex<double> *x_psy;
    complex<double> *psy_x;
    complex<double> *x_psy_psy_x;

    void renew_psy_null_arrays(Hamiltonian &h);
    
    void drop_all_arrays(Hamiltonian &h);
    
    void get_psy(Hamiltonian &h,
                 double time,
                 complex<double> *psy_in,
                 complex<double> *psy_out);
                 
    void print_psy_quadr(Hamiltonian &h);

    double get_norm2_of_psy(Hamiltonian &h,
                            complex<double> *v);
							
    complex<double> get_c_n(int n_c_n, 
                            double u,
                            double eigen_min, 
                            double eigen_max, 
                            double time);

    void get_x_psy_psy_x(Hamiltonian &h);

private:
    IO_data io_data;
};

#define Psy_H
#endif // !Psy_H




