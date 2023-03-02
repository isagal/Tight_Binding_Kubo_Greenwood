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

#ifndef Point_defects_H

class Point_defects
{

public:
    Point_defects();
    ~Point_defects();

    double rs_weak_conc;
    double rs_strong_conc;
    double lrdp_weak_conc;
    double ci_positive_conc;
    double ci_negative_conc;
    double lrdp_strong_conc;
    double oxygen_graph_double_conc;
    double oxygen_graph_double2_conc;
    double oxygen_graph_double3_conc;
    double oxygen_monomeric_conc;
    double oxygen_dimeric_conc;

    double potential;
    double d;
   	
    void get_input_point_defects_data(); 

    void add_point_defects(int n,
                           int n_x,
                           int n_y,
                           complex<double> *h_1,
                           complex<double> *h_2_u,
                           complex<double> *h_2_d,
                           complex<double> *h_3_u,
                           complex<double> *h_3_d,
                           complex<double> *x,
                           complex<double> *y);

    void add_rs_weak(int n,
                     complex<double> *h_1);

    void add_rs_strong(int n,
                       complex<double> *h_1);

    void add_lrdp_weak(int n,
                       complex<double> *h_1,
                       complex<double> *x,
                       complex<double> *y); 

    void add_ci(double ci_sign,
                int n,
                int n_x,
                complex<double> *h_1,
                complex<double> *x,
                complex<double> *y);
 
    void add_lrdp_strong(int n,
                         complex<double> *h_1,
                         complex<double> *x,
                         complex<double> *y); 

    void add_oxygen_graphitic_double(int n,
                                     int n_x,
                                     int n_y,
                                     complex<double> *h_1,
                                     complex<double> *h_2_u,
                                     complex<double> *h_2_d,
                                     complex<double> *h_3_u,
                                     complex<double> *h_3_d,
                                     complex<double> *x,
                                     complex<double> *y); 

    void add_oxygen_graphitic_double2(int n,
                                      int n_x,
                                      int n_y,
                                      complex<double> *h_1,
                                      complex<double> *x,
                                      complex<double> *y); 

    void add_oxygen_graphitic_double3(int n,
                                      int n_x,
                                      int n_y,
                                      complex<double> *h_1,
                                      complex<double> *x,
                                      complex<double> *y); 

    void add_oxygen_monomeric(int n,
                              int n_x,
                              int n_y,
                              complex<double> *h_1,
                              complex<double> *h_2_u,
                              complex<double> *h_2_d,
                              complex<double> *h_3_u,
                              complex<double> *h_3_d,
                              complex<double> *x,
                              complex<double> *y); 

    void add_oxygen_dimeric(int n,
                            int n_x,
                            int n_y,
                            complex<double> *h_1,
                            complex<double> *h_2_u,
                            complex<double> *h_2_d,
                            complex<double> *h_3_u,
                            complex<double> *h_3_d,
                            complex<double> *x,
                            complex<double> *y); 

    double get_ditance(double x1,
                       double x2,
                       double y1,
                       double y2,
                       double cutoff_r); 

    double get_ditance_for_ci(double x1,
                              double x2,
                              double y1,
                              double y2,
                              double height); 

    double get_top_or_down_site(int n,
                                int n_x,
                                int n_y,
                                int n_site);

    double get_random_site_not_close_to_borders(int n,
                                                int n_x,
                                                int n_y);

private:
    IO_data io_data;
};

#define Point_defects_H
#endif // !Point_defects_H




