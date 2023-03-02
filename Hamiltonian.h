#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <complex>

#define MKL_Complex16 std::complex<double>
#define MKL_Complex8  std::complex<float>
#include <mkl.h>
#include "mkl_spblas.h"

#include "IO_data.h"
#include "Strain.h"
#include "Point_defects.h"
#include "Linear_defects.h"


using namespace std;

#ifndef Hamiltonian_H

class Hamiltonian
{

public:
    Hamiltonian();
    ~Hamiltonian();

    int n, n_x, n_y, dn_x, dn_y, dn_x_dos, dn_y_dos;

    double u;
    double emax, emin;
    int ndiag;

    int cont_frac_n, n_energy_points, n_realizations;
    double smooth;

    complex<double> *h_1, *h_2_u, *h_2_d, *h_3_u, *h_3_d;
    complex<double> *x, *y;
    complex<double> *diags, *diags_norm;

    int *idiag;

    complex<double> *alpha, *beta;

    double *energy_array;

    void drop_all_arrays();

    void get_input_lattice_data();
    
    void get_lattice();
    
    void initialize_hamiltonian();
    
    void get_emax_emin();
    
    void normalize_hamiltonian();
    
    void initialize_x_y();
    
    void introduce_impurities();
    
    void introduce_strain();
    
    void get_sparce_matrixes();
    

    double get_distance_between_points(double x1,
                                       double x2,
                                       double y1,
                                       double y2);
                                       
    void mult_matrix(complex<double> *in,
                     complex<double> *out,
                     complex<double> *d);

    void print_h_1();
	
    void get_alpha_beta(complex<double> *starting_array);
	
    double get_continued_fraction (int cfn,
                                   double energy);

    void get_energy_array(double energy_begin,
                          double energy_end);

    complex<double> get_complex_dot (complex<double> *a,
                                     complex<double> *b);
	   				      			 
    void multiply_array_with_scalar (complex<double> *out, 
                                     complex<double> *in, 
                                     complex<double> number);
								     
    void subtract_array_x_form_y (complex<double> *y, 
                                  complex<double> *x);
							      
    void add_array_x_to_y (complex<double> *y, 
                           complex<double> *x);
						   
    Point_defects pd;
    Linear_defects ld;
    Strain strain;

private:
    IO_data io_data;
};

#define Hamiltonian_H
#endif // !Hamiltonian_H




