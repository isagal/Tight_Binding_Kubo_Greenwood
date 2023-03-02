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
#include "Watch.h"
#include "Psy.h"
#include "Hamiltonian.h"


using namespace std;

#ifndef KuboGreenwood_H

class KuboGreenwood
{

public:
    KuboGreenwood();
    ~KuboGreenwood();

    int flag_dos, flag_diff, flag_cond;

    double energy_start, energy_end;

    double diff_start_time, diff_step_time;
    int diff_n_steps;

    double cond_time;

    double *dos;	
    double *el_conc;

    double *diff;	
    double *cond;
    double *cond_in_s;

    void get_input_kb_data();  
    void drop_all_arrays();
    void write_cont_frac_test();
    void calculate_kb();
    void calculate_dos();
    void calculate_diffusivity();
    void calculate_conductivity();
    void print_dos();

    double calculate_dos_for_selected_energy(double energy);
    
    double calculate_diff_for_selected_energy_and_time(double energy,
                                                       double time,
                                                       double &norma);
                                                       
    double calculate_cond_for_selected_energy_and_time(double energy,
                                                       double time,
                                                       double &norma);
                                                       
    void get_psy_propagation(double time, 
                             double &norma);

    void print_diff(int n_time_step,
                    double time);
                    
    void print_cond();

    void calculate_el_conc();

    void get_min_cond_in_s();

private:
    IO_data io_data; Hamiltonian ham; Psy psy; Watch watch;
};

#define KuboGreenwood_H
#endif // !KuboGreenwood_H




