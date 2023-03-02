#include "Hamiltonian.h"

Hamiltonian::Hamiltonian() : io_data(), pd(), ld(), strain()
{
    u = 2.78;
	
    ndiag = 5;    //Number of diagonals

    get_input_lattice_data();

    dn_x_dos = int(n_x * 0.9);
    dn_y_dos = int(n_y * 0.9);
	
    x = y = NULL;
    h_1 = h_2_u = h_2_d = h_3_u = h_3_d = NULL;	
    diags = diags_norm = NULL;
    idiag = NULL;

    h_1   = new complex<double>[n];
    h_2_u = new complex<double>[n];
    h_2_d = new complex<double>[n];
    h_3_u = new complex<double>[n];
    h_3_d = new complex<double>[n];

    x = new complex<double>[n]; 
    y = new complex<double>[n]; 

    diags_norm = new complex<double>[ndiag*n];
    diags = new complex<double>[ndiag*n];
    idiag = new int[ndiag];

    alpha = new complex<double>[cont_frac_n];
    beta  = new complex<double>[cont_frac_n];

    energy_array = new double[n_energy_points];
}

/*----------------------------------------------------------------------------*/

Hamiltonian::~Hamiltonian()
{
    if (x != NULL) {delete [] x;}
    if (y != NULL) {delete [] y;}

    if (h_1   != NULL) {delete [] h_1  ;}
    if (h_2_u != NULL) {delete [] h_2_u;}
    if (h_2_d != NULL) {delete [] h_2_d;}
    if (h_3_u != NULL) {delete [] h_3_u;}
    if (h_3_d != NULL) {delete [] h_3_d;}

    if (idiag      != NULL) {delete [] idiag     ;}
    if (diags      != NULL) {delete [] diags     ;}
    if (diags_norm != NULL) {delete [] diags_norm;}

    if (alpha != NULL) {delete [] alpha;}
    if (beta  != NULL) {delete [] beta ;}

    if (energy_array != NULL) {delete [] energy_array;}
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::get_input_lattice_data() 
{
    //gets all data form the lattice input file 
    string s;
    fstream input_file("./input/input_lattice.txt");

    io_data.write_line_to_output_file("Input lattice data:\n");

    getline(input_file, s);
    n_x = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Number of lattice sites along X :") + to_string(n_x));

    getline(input_file, s);
    n_y = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Number of lattice sites along Y :") + to_string(n_y));

    n = n_x * n_y;
    io_data.write_line_to_output_file(string("Total number of lattice sites   :") + to_string(n));

    getline(input_file, s);
    dn_x = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Number of wave-packet sites along X :") + to_string(dn_x));

    getline(input_file, s);
    dn_y = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Number of wave-packet sites along Y :") + to_string(dn_y));

    getline(input_file, s);
    cont_frac_n = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Continued fraction truncation term :") + to_string(cont_frac_n));

    getline(input_file, s);
    smooth = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Smoothing factor for continued fraction :") + to_string(smooth));

    getline(input_file, s);
    n_energy_points = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Number of energy points :") + to_string(n_energy_points));

    getline(input_file, s);
    n_realizations = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Number of realizations :") + to_string(n_realizations));

    io_data.write_line_to_output_file("____________________\n"); 

    input_file.close();
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::drop_all_arrays()
{
    bzero(x,     sizeof(complex<double>) * n);
    bzero(y,     sizeof(complex<double>) * n);
    bzero(h_1,   sizeof(complex<double>) * n);
    bzero(h_2_u, sizeof(complex<double>) * n);
    bzero(h_2_d, sizeof(complex<double>) * n);
    bzero(h_3_u, sizeof(complex<double>) * n);
    bzero(h_3_d, sizeof(complex<double>) * n);

    bzero(diags     , sizeof(complex<double>) * ndiag * n);
    bzero(diags_norm, sizeof(complex<double>) * ndiag * n);

    bzero(idiag, sizeof(int) * ndiag);

    bzero(alpha, sizeof(complex<double>) * cont_frac_n);
    bzero(beta,  sizeof(complex<double>) * cont_frac_n);

    bzero(energy_array, sizeof(double) * n_energy_points);
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::get_lattice()
{
    initialize_hamiltonian();
    initialize_x_y();
    introduce_strain();
    introduce_impurities();
    get_emax_emin();
    get_sparce_matrixes();
    normalize_hamiltonian();
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::initialize_hamiltonian()
{
    idiag[0] = -n_x; 
    idiag[1] = -1; 
    idiag[2] =  0; 
    idiag[3] =  1; 
    idiag[4] =  n_x;
        
    double gamma = 1.0;	//hopping integral in |u| units

    int coin = 0;

    for(int i = 0 ; i < (n - 1); i++)    //initializing of the h_2_u
    { 
        if ((i + 1) % n_x != 0) {h_2_u[i] = -gamma;}
    }  
                  
    for(int i = 0 ; i < n - n_x ; i++)    //initializing of the h_3_u
    {
        if ((coin % 2 == 0) && ((i + 1) % 2 != 0)) {h_3_u[i] = -gamma;}
        if ((coin % 2 != 0) && (i % 2 != 0))       {h_3_u[i] = -gamma;}                      
        if ((i + 1) % n_x == 0)                    {coin += 1;}     
    }

    for(int i = 1   ; i < n ; i++) {h_2_d[i] = h_2_u[i - 1  ];}    //initializing of the h_2_d
    for(int i = n_x ; i < n ; i++) {h_3_d[i] = h_3_u[i - n_x];}    //initializing of the h_3_d
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::get_sparce_matrixes()
{
    memcpy(&diags[0 * n], h_3_d, sizeof(complex<double>) * n);
    memcpy(&diags[1 * n], h_2_d, sizeof(complex<double>) * n);
    memcpy(&diags[2 * n], h_1  , sizeof(complex<double>) * n);
    memcpy(&diags[3 * n], h_2_u, sizeof(complex<double>) * n);
    memcpy(&diags[4 * n], h_3_u, sizeof(complex<double>) * n);

    memcpy(&diags_norm[0 * n], h_3_d, sizeof(complex<double>) * n);
    memcpy(&diags_norm[1 * n], h_2_d, sizeof(complex<double>) * n);
    memcpy(&diags_norm[2 * n], h_1  , sizeof(complex<double>) * n);
    memcpy(&diags_norm[3 * n], h_2_u, sizeof(complex<double>) * n);
    memcpy(&diags_norm[4 * n], h_3_u, sizeof(complex<double>) * n);
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::get_emax_emin()
{
    double *h_1_double = new double[n];
    double *h_2_u_abs  = new double[n];
    double *h_2_d_abs  = new double[n];
    double *h_3_u_abs  = new double[n];
    double *h_3_d_abs  = new double[n];
    double *all_h_abs  = new double[n];
    double *temp_array = new double[n];

    bzero(h_1_double, sizeof(double) * n);
    bzero(h_2_u_abs,  sizeof(double) * n);
    bzero(h_2_d_abs,  sizeof(double) * n);
    bzero(h_3_u_abs,  sizeof(double) * n);
    bzero(h_3_d_abs,  sizeof(double) * n);
    bzero(all_h_abs,  sizeof(double) * n);
    bzero(temp_array, sizeof(double) * n);

    vzAbs(n, h_3_d, h_3_d_abs);
    vzAbs(n, h_2_d, h_2_d_abs);
    vzAbs(n, h_2_u, h_2_u_abs);
    vzAbs(n, h_3_u, h_3_u_abs);

    for(int i = 0 ; i < n ; i++) {h_1_double[i] = real(h_1[i]);}

    vdAdd(n, h_3_d_abs, h_3_u_abs, all_h_abs);
    vdAdd(n, all_h_abs, h_2_u_abs, all_h_abs);
    vdAdd(n, all_h_abs, h_2_d_abs, all_h_abs);

    vdAdd(n, h_1_double, all_h_abs, temp_array);

    emax = temp_array[0];
    for(int i = 1 ; i < n ; i++)
    {
        if(emax < temp_array[i]) {emax = temp_array[i];}
    }

    vdSub(n, h_1_double, all_h_abs, temp_array);
    emin = temp_array[0];
    for(int i = 1 ; i < n ; i++)
    {
        if(emin > temp_array[i]) {emin = temp_array[i];}
    }

    emin *= 1.1; 
    emax *= 1.1; 

    if (h_1_double != NULL) delete [] h_1_double;
    if (h_2_u_abs  != NULL) delete [] h_2_u_abs;
    if (h_2_d_abs  != NULL) delete [] h_2_d_abs;
    if (h_3_u_abs  != NULL) delete [] h_3_u_abs;
    if (h_3_d_abs  != NULL) delete [] h_3_d_abs;
    if (all_h_abs  != NULL) delete [] all_h_abs;
    if (temp_array != NULL) delete [] temp_array;
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::normalize_hamiltonian()
{
    double eigen_sum  = emax + emin;
    double eigen_diff = emax - emin;

    for(int i = 0 * n ; i < 2 * n ; i++) {diags_norm[i] =  2.0 * diags_norm[i] / eigen_diff;}
    for(int i = 2 * n ; i < 3 * n ; i++) {diags_norm[i] = (2.0 * diags_norm[i] - eigen_sum) / eigen_diff;}
    for(int i = 3 * n ; i < 5 * n ; i++) {diags_norm[i] =  2.0 * diags_norm[i] / eigen_diff;}
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::initialize_x_y()
{
    int ii = 0;
    for(int i = 0 ; i < n ; i++)
    {
        x[i] = ii * sqrt(3.0) / 2.0;
        if(ii == (n_x - 1)) {ii  = 0;}
        else                {ii += 1;}
    }    

    ii = 0;
    int index_up_down = 0;    //up down index for graphene lattice, temporary variable
    for(int i = 0 ; i < n ; i++)
    {
        if(i % n_x == 0) 
        {
			ii += 1;
            index_up_down += 1;
        }
        y[i] = (ii - 1.0) * 1.5;
        if(index_up_down % 2 != 0) 
        {
            if(i % 2 == 0) {y[i] += 0.5;}
        }
        else
        {
            if(i % 2 != 0) {y[i] += 0.5;}
        }
    }	
}

/*----------------------------------------------------------------------------*/

double Hamiltonian::get_distance_between_points(double x1,
												double x2,
												double y1,
												double y2)
{
    double dx = x1 - x2;
    double dy = y1 - y2;
    return sqrt(dx*dx + dy*dy);
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::introduce_impurities()
{
    pd.add_point_defects(n, n_x, n_y, h_1, h_2_u, h_2_d, h_3_u, h_3_d, x, y);
    ld.add_linear_defects(n, n_x, h_1, x, y);
    print_h_1();
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::introduce_strain()
{
    strain.add_strain(n, n_x, h_2_u, h_2_d, h_3_u, h_3_d, x, y);
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::print_h_1()
{
    double re = 0.0;
    ofstream fres("./output/h_1.dat", ios_base::out);

    for(int i = 0 ; i < n ; i++)
    {
        re = real(h_1[i]);
        if(re < 1e-10) {re = 0.0;}	
        if((i + 1) % n_x == 0) {fres << real(h_1[i]) << endl;}
        else {fres << real(h_1[i]) << " ";}
    }	
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::get_alpha_beta (complex<double> *starting_array)
{
    /**
     *          a1  b1  0  0  0  0  0  ....      
     *          b1  a2  b2 0  0  0  0  ....     
     *          0   b2  a3 b3 0  0  0  ....    
     * H_try=   0   0   b3 a4 b4 0  0  ....     
     *          ...........................        
     *          .....................b(n-1)       
     *          0  0  0  .......  b(n-1) an   
     */

    complex<double> *gamma        = new complex<double> [n];
    complex<double> *gamma_old    = new complex<double> [n];
    complex<double> *gamma_new    = new complex<double> [n];
    complex<double> *h_mult_gamma = new complex<double> [n];

    complex<double> *gamma_mult_alpha    = new complex<double> [n];
    complex<double> *gamma_old_mult_beta = new complex<double> [n];

    bzero(gamma,               sizeof(complex<double>) * n);
    bzero(gamma_old,           sizeof(complex<double>) * n);
    bzero(gamma_new,           sizeof(complex<double>) * n);
    bzero(h_mult_gamma,        sizeof(complex<double>) * n);
    bzero(gamma_mult_alpha,    sizeof(complex<double>) * n);
    bzero(gamma_old_mult_beta, sizeof(complex<double>) * n);

    memcpy(gamma, starting_array, sizeof(complex<double>) * n);

    bzero(alpha, sizeof(complex<double>) * cont_frac_n);
    bzero(beta,  sizeof(complex<double>) * cont_frac_n);

    for(int i = 0 ; i < cont_frac_n ; i++)
    {
        mult_matrix (gamma, h_mult_gamma, diags);
        alpha[i] += get_complex_dot(gamma, h_mult_gamma);
        multiply_array_with_scalar(gamma_mult_alpha, gamma, alpha[i]);

        bzero(gamma_new, sizeof(complex<double>) * n);
        add_array_x_to_y (gamma_new, h_mult_gamma);
        subtract_array_x_form_y	(gamma_new, gamma_mult_alpha);	
        if(i > 0)
        {
            multiply_array_with_scalar (gamma_old_mult_beta, gamma_old, beta[i - 1]);
            subtract_array_x_form_y	(gamma_new, gamma_old_mult_beta);					
        }      

        beta[i] += beta[i] + get_complex_dot(gamma_new, gamma_new);
        beta[i]  = sqrt(beta[i]);
        memcpy(gamma_old, gamma, sizeof(complex<double>) * n);
		
        multiply_array_with_scalar (gamma, gamma_new, (1.0 / beta[i]));	
    }

    if (gamma        != NULL) delete [] gamma;
    if (gamma_old    != NULL) delete [] gamma_old;
    if (gamma_new    != NULL) delete [] gamma_new;
    if (h_mult_gamma != NULL) delete [] h_mult_gamma;

    if (gamma_mult_alpha    != NULL) delete [] gamma_mult_alpha;
    if (gamma_old_mult_beta != NULL) delete [] gamma_old_mult_beta;
}

/*----------------------------------------------------------------------------*/

double Hamiltonian::get_continued_fraction (int cfn, double energy)
{
    complex<double> k1 (0.0, 0.0);               
    complex<double> kn (0.0, 0.0);  
    complex<double> e1 (energy,smooth);                 //temporary variable
    complex<double> e2 = e1 - alpha[cfn - 1];   //temporary variable
    complex<double> last_beta_sqr = beta[cfn-1] * beta[cfn-1];

    complex<double> temp_var = sqrt(complex<double>(4.0, 0.0) * last_beta_sqr - e2 * e2);
    complex<double> sum_e(0.0, 0.0);

    sum_e  = e2 - complex<double>(0.0,1.0) * temp_var;
    sum_e /= complex<double> (2.0, 0.0) * last_beta_sqr;

    kn = e2 - last_beta_sqr * sum_e;

    for(int i = (cfn - 2) ; i > -1 ; i -= 1)
    {
        k1 = e1 - alpha[i] - beta[i] * beta[i] / kn;
        kn = k1;
    }

    return imag(complex<double>(1.0,0.0) / k1);
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::get_energy_array(double energy_begin,
                                   double energy_end)
{
    for(int i = 0 ; i < n_energy_points ; i++)
    {
        energy_array[i]  = energy_begin + i * energy_end / ((n_energy_points - 1.0) / 2.0);
    }
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::mult_matrix (complex<double> *in,
                               complex<double> *out,
                               complex<double> *d)
{
    //d - diags matrix
    char transa= 'N';
    mkl_zdiagemv(&transa, &n, d, &n, idiag, &ndiag, in, out);
}

/*----------------------------------------------------------------------------*/

complex<double> Hamiltonian::get_complex_dot (complex<double> *a,
                                              complex<double> *b)
{
    complex<double> result=0.0;
    int incx = 1, incy = 1;
    zdotc(&result, &n, a, &incx, b, &incy);
    return result;
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::multiply_array_with_scalar (complex<double> *out, 
                                              complex<double> *in, 
                                              complex<double> number)
{
    complex<double> complex_zero (0.0, 0.0);
    cblas_zaxpby (n, &number, in, 1, &complex_zero, out, 1);
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::subtract_array_x_form_y (complex<double> *y, 
                                           complex<double> *x)
{
    // y = (-1) * x + y
    complex<double> minus_one (-1.0, 0.0); 
    cblas_zaxpy (n, &minus_one, x, 1, y, 1);
}

/*----------------------------------------------------------------------------*/

void Hamiltonian::add_array_x_to_y (complex<double> *y, 
                                    complex<double> *x)
{
    // y = x + y
    complex<double> one (1.0, 0.0); 
    cblas_zaxpy (n, &one, x, 1, y, 1);
}







