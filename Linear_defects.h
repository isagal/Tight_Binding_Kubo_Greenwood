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

#ifndef Linear_defects_H

class Linear_defects
{

public:
    Linear_defects();
    ~Linear_defects();

    int line_total_number, line_potential_sign, line_orientation;
    double line_corr_angle;
	
    void get_input_linear_defects_data();
    
    void add_linear_defects(int n,
                            int n_x,
                            complex<double> *h_1,
                            complex<double> *x,
                            complex<double> *y);
                            
    void add_strictly_zigzag_lines(int n,
                                   complex<double> *h_1,
                                   complex<double> *x,
                                   complex<double> *y);
							  	   
    void add_strictly_armchair_lines(int n,
                                     complex<double> *h_1,
                                     complex<double> *x,
                                     complex<double> *y); 
							  		 
    void add_random_lines(int n,
                          complex<double> *h_1,
                          complex<double> *x,
                          complex<double> *y); 
						  
    void add_correlated_lines(int n,
                              complex<double> *h_1,
                              complex<double> *x,
                              complex<double> *y); 
							  

    double get_ditance(double x1,
                       double x2,
                       double y1,
                       double y2); 

    double get_potential_for_line(); 
	
    double calculate_potential_from_line_to_site(double r,
                                                 double uu);
					   							 
    void get_random_point_for_future_line(int n, 
                                          complex<double> *x,
                                          complex<double> *y,
                                          double &x1, 
                                          double &y1);
                                          
    void get_k_and_b_for_a_random_line(int n, 
                                       complex<double> *x,
                                       complex<double> *y,
                                       double &k, 
                                       double &b);
                                       
    double get_distance_from_line_to_site(complex<double> *x,
                                          complex<double> *y,
                                          double k,
                                          double b,
                                          int site_number);
                                          
    void get_k_and_b_for_a_correlated_line(int n,
                                           complex<double> *x,
                                           complex<double> *y,
                                           double k_first,
                                           double b_first,
                                           double &k,
                                           double &b);
                                           
    void add_first_line_correlated_case(int n,
                                        complex<double> *x,
                                        complex<double> *y,
                                        double &k_first,
                                        double &b_first); 

private:
    IO_data io_data;
};

#define Linear_defects_H
#endif // !Linear_defects_H




