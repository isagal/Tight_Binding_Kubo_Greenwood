#include "Strain.h"

Strain::Strain() : io_data()
{
    /**
     * Bond vectors followed below:
     *     
     *      delta_2
     *         |
     *        / \
     * delta_3   delta_1
     * 
     */

    get_input_strain_data(); 
    //Poisson ratio is 0.15 (between ones for graphite 0.165 and graphene 0.14)
    poisson = 0.15;
}

/*----------------------------------------------------------------------------*/

Strain::~Strain()
{

}

/*----------------------------------------------------------------------------*/

void Strain::get_input_strain_data() 
{
    //gets all data form the point defects input file 
    string s;
    fstream input_file("./input/input_strain.txt");

    io_data.write_line_to_output_file("Input strain data:\n");

    getline(input_file, s);
    uniax_flag = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file("Uniaxial strain (0 - no, 1 - x, 2 - y)      : " + to_string(uniax_flag));

    getline(input_file, s);
    uniax_val = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Uniaxial strain value (from 0.0 to 1.0)     : " + to_string(uniax_val));

    getline(input_file, s);
    shear_flag = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file("Shear strain (0 - no, 1 - x, 2 - y, 3 - xy) : " + to_string(shear_flag));

    getline(input_file, s);
    shear_val = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Shear strain value (from 0.0 to 1.0)        : " + to_string(shear_val));

    io_data.write_line_to_output_file("____________________\n"); 

    input_file.close();
}

/*----------------------------------------------------------------------------*/

void Strain::add_strain(int n,
                        int n_x,
                        complex<double> *h_2_u,
                        complex<double> *h_2_d,
                        complex<double> *h_3_u,
                        complex<double> *h_3_d,
                        complex<double> *x,
                        complex<double> *y) 
{
    if(uniax_flag != 0)
    {
        uniax_modify_hamiltonian(n,
                                 h_2_u,
                                 h_2_d,
                                 h_3_u,
                                 h_3_d); 
        uniax_modify_x_y(n,
                         x,
                         y);
    }
    if(shear_flag != 0)
    {
        shear_modify_hamiltonian(n,
                                 n_x,
                                 h_2_u,
                                 h_2_d,
                                 h_3_u,
                                 h_3_d,
                                 x,
                                 y); 
        shear_modify_x_y(n,
                         x,
                         y);
    }
}

/*----------------------------------------------------------------------------*/

void Strain::uniax_modify_hamiltonian(int n,
                                      complex<double> *h_2_u,
                                      complex<double> *h_2_d,
                                      complex<double> *h_3_u,
                                      complex<double> *h_3_d) 
{
    double p = poisson; 
    double d = uniax_val;

    double gamma_0_delta_2   = 0.0;
    double gamma_0_delta_1_3 = 0.0;
        
    if(uniax_flag == 1)
    {
        gamma_0_delta_1_3 = sqrt((1.0 + d) * (1.0 + d) * 3.0 / 4.0 + (1.0 - d * p) * (1.0 - d * p) * 1.0 / 4.0);
        gamma_0_delta_2   = 1.0 - d * p;          
    }
    if(uniax_flag == 2)
    {
        gamma_0_delta_1_3 = sqrt((1.0 - d * p) * (1.0 - d * p) * 3.0 / 4.0 + (1.0 + d) * (1.0 + d) * 1.0 / 4.0);
        gamma_0_delta_2   = 1.0 + d;
    }

    complex<double> scale_coef_2   (exp(-3.37 * (gamma_0_delta_2   - 1.0)),0.0);
    complex<double> scale_coef_1_3 (exp(-3.37 * (gamma_0_delta_1_3 - 1.0)),0.0);

    cblas_zscal(n, &scale_coef_2,   h_3_u, 1);
    cblas_zscal(n, &scale_coef_2,   h_3_d, 1);
    cblas_zscal(n, &scale_coef_1_3, h_2_u, 1);
    cblas_zscal(n, &scale_coef_1_3, h_2_d, 1);
}

/*----------------------------------------------------------------------------*/

void Strain::uniax_modify_x_y(int n,
                              complex<double> *x,
                              complex<double> *y) 
{
    double p = poisson; 
    double d = uniax_val;

    complex<double> scale_coef_1 ((1.0 + d)    ,0.0);
    complex<double> scale_coef_2 ((1.0 - d * p),0.0);

    if(uniax_flag == 1)	
    {
        cblas_zscal(n, &scale_coef_1, x, 1);
        cblas_zscal(n, &scale_coef_2, y, 1);
    }
    if(uniax_flag == 2)
    {
        cblas_zscal(n, &scale_coef_2, x, 1);
        cblas_zscal(n, &scale_coef_1, y, 1);
    }
}

/*----------------------------------------------------------------------------*/

void Strain::shear_modify_hamiltonian(int n,
                                      int n_x,
                                      complex<double> *h_2_u,
                                      complex<double> *h_2_d,
                                      complex<double> *h_3_u,
                                      complex<double> *h_3_d,
                                      complex<double> *x,
                                      complex<double> *y) 
{
    double p = poisson; 
    double d = shear_val;

    double new_x_0   = 0.0;
    double new_y_0   = 0.0;
    double new_x_1   = 0.0;
    double new_y_1   = 0.0;
    double new_x_2   = 0.0;
    double new_y_2   = 0.0;
    double new_x_n_x = 0.0;
    double new_y_n_x = 0.0;

    double delta_1_old = get_distance_between_two_points(real(x[1]),   real(y[1]),   real(x[0]), real(y[0]));
    double delta_2_old = get_distance_between_two_points(real(x[n_x]), real(y[n_x]), real(x[0]), real(y[0]));
    double delta_3_old = get_distance_between_two_points(real(x[2]),   real(y[2]),   real(x[1]), real(y[1]));
                
    if(shear_flag == 1)
    {
        new_x_0   = real(x[0])   + d * real(y[0]);
        new_y_0   = real(y[0]);
        new_x_1   = real(x[1])   + d * real(y[1]);
        new_y_1   = real(y[1]);
        new_x_2   = real(x[2])   + d * real(y[2]);
        new_y_2   = real(y[2]);
        new_x_n_x = real(x[n_x]) + d * real(y[n_x]);
        new_y_n_x = real(y[n_x]);
    }
    if(shear_flag == 2)
    {
        new_x_0   = real(x[0]);
        new_y_0   = real(y[0])   + d * real(x[0]);
        new_x_1   = real(x[1]);
        new_y_1   = real(y[1])   + d * real(x[1]);
        new_x_2   = real(x[2]);
        new_y_2   = real(y[2])   + d * real(x[2]);
        new_x_n_x = real(x[n_x]);
        new_y_n_x = real(y[n_x]) + d * real(x[n_x]);
    }
    if(shear_flag == 3)
    {
        new_x_0   = real(x[0])   + d * real(y[0]);
        new_y_0   = real(y[0])   + d * real(x[0]);
        new_x_1   = real(x[1])   + d * real(y[1]);
        new_y_1   = real(y[1])   + d * real(x[1]);
        new_x_2   = real(x[2])   + d * real(y[2]);
        new_y_2   = real(y[2])   + d * real(x[2]);
        new_x_n_x = real(x[n_x]) + d * real(y[n_x]);
        new_y_n_x = real(y[n_x]) + d * real(x[n_x]);
    }
        
    double delta_1_new = get_distance_between_two_points(new_x_1,   new_y_1,   new_x_0, new_y_0);
    double delta_2_new = get_distance_between_two_points(new_x_n_x, new_y_n_x, new_x_0, new_y_0);
    double delta_3_new = get_distance_between_two_points(new_x_2,   new_y_2,   new_x_1, new_y_1);
        
    double gamma_0_delta_1 = delta_1_new / delta_1_old;
    double gamma_0_delta_2 = delta_2_new / delta_2_old;
    double gamma_0_delta_3 = delta_3_new / delta_3_old;

    int coin = 0;

    complex<double> scale_coef_1 (exp(-3.37 * (gamma_0_delta_1 - 1.0)),0.0);
    complex<double> scale_coef_2 (exp(-3.37 * (gamma_0_delta_2 - 1.0)),0.0);
    complex<double> scale_coef_3 (exp(-3.37 * (gamma_0_delta_3 - 1.0)),0.0);

    cblas_zscal(n, &scale_coef_2,   h_3_u, 1);
    cblas_zscal(n, &scale_coef_2,   h_3_d, 1);

    for(int i = 0 ; i < (n - 1) ; i++)
    {
        if(coin % 2 == 0)
        {
            if(i % 2 == 0) {h_2_u[i] *= scale_coef_1;} 
            else           {h_2_u[i] *= scale_coef_3;}
        }
        else
        {
            if(i % 2 == 0) {h_2_u[i] *= scale_coef_3;} 
            else           {h_2_u[i] *= scale_coef_1;}
        }           
        if((i + 1) % n_x == 0) {coin += 1;}	
    }
    memcpy(&h_2_d[1], h_2_u, sizeof(complex<double>) * (n - 1));
}

/*----------------------------------------------------------------------------*/

void Strain::shear_modify_x_y(int n,
                              complex<double> *x,
                              complex<double> *y) 
{
    double p = poisson; 
    double d = shear_val;

    if(shear_flag == 1)	
    {
        for(int i = 0 ; i < n ; i++) {x[i] += y[i] * d;}
    }
    if(shear_flag == 2)
    {
        for(int i = 0 ; i < n ; i++) {y[i] += x[i] * d;}
    }
    if(shear_flag == 3)
    {
        complex<double> *old_x = new complex<double> [n];
        complex<double> *old_y = new complex<double> [n];
        bzero(old_x, sizeof(complex<double>) * n);
        bzero(old_y, sizeof(complex<double>) * n);
        memcpy(old_x, x, sizeof(complex<double>) * n);		
        memcpy(old_y, y, sizeof(complex<double>) * n);	

        for(int i = 0 ; i < n ; i++) 
        {
            x[i] += old_y[i] * d;
            y[i] += old_x[i] * d;			
        }
	
        if (old_x != NULL) delete [] old_x;
        if (old_y != NULL) delete [] old_y;
    }
}

/*----------------------------------------------------------------------------*/

double Strain::get_distance_between_two_points(double x1,
                                               double y1,
                                               double x2,
                                               double y2) 
{
    double xx = (x2 - x1) * (x2 - x1);
    double yy = (y2 - y1) * (y2 - y1);
    return sqrt(xx + yy);
}







