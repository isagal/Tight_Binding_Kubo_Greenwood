#include "Linear_defects.h"

Linear_defects::Linear_defects() : io_data()
{
    get_input_linear_defects_data(); 
}

/*----------------------------------------------------------------------------*/

Linear_defects::~Linear_defects()
{

}

/*----------------------------------------------------------------------------*/

void Linear_defects::get_input_linear_defects_data() 
{
    //gets all data form the point defects input file 
    string s;
    fstream input_file("./input/input_linear_defects.txt");

    io_data.write_line_to_output_file("Input linear defects data:\n");

    getline(input_file, s);
    line_total_number = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file("Number of lines :" + to_string(line_total_number));

    getline(input_file, s);
    line_potential_sign = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file("Type of the line potential :" + to_string(line_potential_sign) );

    getline(input_file, s);
    line_orientation = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file("Orientation of the first line :" + to_string(line_orientation));

    getline(input_file, s);
    line_corr_angle = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Correlation angle in Pi units :" + to_string(line_corr_angle));

    io_data.write_line_to_output_file("____________________\n"); 

    input_file.close();
}

/*----------------------------------------------------------------------------*/

void Linear_defects::add_linear_defects(int n,
                                        int n_x,
                                        complex<double> *h_1,
                                        complex<double> *x,
                                        complex<double> *y) 
{
    if(line_total_number > 0)
    {
        if(line_orientation == 0) 
        {
            add_strictly_zigzag_lines(n, h_1, x, y);
        }
        if(line_orientation == 1)
        {
            add_strictly_armchair_lines(n, h_1, x, y);
        }
        if(line_orientation == 2)
        {
            add_random_lines(n, h_1, x, y);
        }
        if((line_orientation == 3) || (line_orientation == 4) || (line_orientation == 5))
        {
            add_correlated_lines(n, h_1, x, y);
        }
    }
}

/*----------------------------------------------------------------------------*/

void Linear_defects::add_strictly_zigzag_lines(int n,
                                               complex<double> *h_1,
                                               complex<double> *x,
                                               complex<double> *y) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    double r = 0.0;
    double potential = 0.0;
    double *y_of_lines = new double [line_total_number];
    double *u_of_lines = new double [line_total_number];

    bzero(y_of_lines, sizeof(double) * line_total_number);
    bzero(u_of_lines, sizeof(double) * line_total_number);

    for(int i = 0 ; i < line_total_number ; i++)
    {
        y_of_lines[i] = (real(y[n - 1]) - real(y[0])) * drand(gen);
        u_of_lines[i] = get_potential_for_line();        
        for(int j = 0 ; j < n ; j++)
        {
            r = fabs(real(y[j]) - y_of_lines[i]);
            potential = calculate_potential_from_line_to_site(r, u_of_lines[i]);
            h_1[j] += complex<double> (potential,0.0);
        }
    }

    if (y_of_lines != NULL) delete [] y_of_lines;
    if (u_of_lines != NULL) delete [] u_of_lines;
}

/*----------------------------------------------------------------------------*/

void Linear_defects::add_strictly_armchair_lines(int n,
                                                 complex<double> *h_1,
                                                 complex<double> *x,
                                                 complex<double> *y) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    double r = 0.0;
    double potential = 0.0;
    double *x_of_lines = new double [line_total_number];
    double *u_of_lines = new double [line_total_number];

    bzero(x_of_lines, sizeof(double) * line_total_number);
    bzero(u_of_lines, sizeof(double) * line_total_number);

    for(int i = 0 ; i < line_total_number ; i++)
    {
        x_of_lines[i] = (real(x[n - 1]) - real(x[0])) * drand(gen);
        u_of_lines[i] = get_potential_for_line();        

        for(int j = 0 ; j < n ; j++)
        {
            r = fabs(real(x[j]) - x_of_lines[i]);
            potential = calculate_potential_from_line_to_site(r, u_of_lines[i]);
            h_1[j] += complex<double> (potential,0.0);
        }
    }

    if (x_of_lines != NULL) delete [] x_of_lines;
    if (u_of_lines != NULL) delete [] u_of_lines;
}

/*----------------------------------------------------------------------------*/

void Linear_defects::add_random_lines(int n,
                                      complex<double> *h_1,
                                      complex<double> *x,
                                      complex<double> *y) 
{
    //All lines are randomly distributed
    double r = 0;
    double *k          = new double [line_total_number];
    double *b          = new double [line_total_number];
    double *u_of_lines = new double [line_total_number];
    bzero(k,          sizeof(double) * line_total_number);
    bzero(b,          sizeof(double) * line_total_number);
    bzero(u_of_lines, sizeof(double) * line_total_number);

    for(int i = 0 ; i < line_total_number ; i++)
    {
        get_k_and_b_for_a_random_line(n, x, y, k[i], b[i]);
        u_of_lines[i] = get_potential_for_line(); 
        for(int j = 0 ; j < n ; j++)
        {
            r = get_distance_from_line_to_site(x, y, k[i], b[i], j);
            h_1[j] += calculate_potential_from_line_to_site(r, u_of_lines[i]);   
        }
    }

    if (k          != NULL) delete [] k;
    if (b          != NULL) delete [] b;
    if (u_of_lines != NULL) delete [] u_of_lines;
}

/*----------------------------------------------------------------------------*/

void Linear_defects::add_correlated_lines(int n,
                                          complex<double> *h_1,
                                          complex<double> *x,
                                          complex<double> *y) 
{
    double r = 0.0;
    double *k          = new double [line_total_number];
    double *b          = new double [line_total_number];
    double *u_of_lines = new double [line_total_number];
    bzero(k,          sizeof(double) * line_total_number);
    bzero(b,          sizeof(double) * line_total_number);
    bzero(u_of_lines, sizeof(double) * line_total_number);

    add_first_line_correlated_case(n, x, y, k[0], b[0]);
            
    for(int j = 0 ; j < line_total_number ; j++)
    {
        get_k_and_b_for_a_correlated_line(n, x, y, k[0], b[0], k[j], b[j]);
        u_of_lines[j] = get_potential_for_line();     
        for(int f = 0 ; f < n ; f++)
        {
            r = get_distance_from_line_to_site(x, y, k[j], b[j], f);
            h_1[f] += calculate_potential_from_line_to_site(r, u_of_lines[j]);
        } 
    }

    if (k          != NULL) delete [] k;
    if (b          != NULL) delete [] b;
    if (u_of_lines != NULL) delete [] u_of_lines;
}

/*----------------------------------------------------------------------------*/

double Linear_defects::get_distance_from_line_to_site(complex<double> *x,
                                                      complex<double> *y,
                                                      double k,
                                                      double b,
                                                      int site_number) 
{
    double k0 = -1.0 / tan(atan(k));
    double b0 = real(y[site_number]) - k0 * real(x[site_number]);
    double x0 = (b0 - b) / (k - k0);
    double y0 = k0 * (b0 - b) / (k - k0) + b0;
    double dx = real(x[site_number]) - x0;
    double dy = real(y[site_number]) - y0;

    return sqrt(dx * dx + dy * dy);
}

/*----------------------------------------------------------------------------*/

void Linear_defects::get_k_and_b_for_a_random_line(int n, 
                                                   complex<double> *x,
                                                   complex<double> *y,
                                                   double &k, 
                                                   double &b)
{
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    get_random_point_for_future_line(n, x, y, x1, y1);
    get_random_point_for_future_line(n, x, y, x2, y2);
    k = (y2 - y1) / (x2 - x1);
    b = y1 - k * x1;
}

/*----------------------------------------------------------------------------*/

void Linear_defects::get_random_point_for_future_line(int n, 
                                                      complex<double> *x,
                                                      complex<double> *y,
                                                      double &x1, 
                                                      double &y1)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    int a  = int((n - 1) * drand(gen));
    x1 = real(x[a]) + 0.000000001 * (drand(gen) - 0.5);
    y1 = real(y[a]) + 0.000000001 * (drand(gen) - 0.5);
}

/*----------------------------------------------------------------------------*/

double Linear_defects::calculate_potential_from_line_to_site(double r,
                                                             double uu) 
{
    return (uu * 1.544 / (0.780 + 0.046 * r * r));
}

/*----------------------------------------------------------------------------*/

double Linear_defects::get_potential_for_line() 
{
    double potential = 0.0;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    if(line_potential_sign == 0) {potential = 0.25 * drand(gen)       ;}
    if(line_potential_sign == 1) {potential = 0.5  * drand(gen) - 0.25;}

    return potential;
}

/*----------------------------------------------------------------------------*/

void Linear_defects::get_k_and_b_for_a_correlated_line(int n,
                                                       complex<double> *x,
                                                       complex<double> *y,
                                                       double k_first,
                                                       double b_first,
                                                       double &k,
                                                       double &b) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    double x1 = 0.0;
    double y1 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    double a  = 0.0;

    get_random_point_for_future_line(n, x, y, x1, y1);
    get_random_point_for_future_line(n, x, y, x2, y2);

    if(line_orientation == 3)
    {
        a = (drand(gen) - 0.5) * line_corr_angle * M_PI;
        k = tan(atan(k_first) + a);
        b = y1 - k * x1;
    }
    if(line_orientation == 4)
    {
        a = (drand(gen) - 0.5) * line_corr_angle * M_PI;
        k = tan(atan(k_first) + a);
        b = y1 - k * x1;
    }
    if(line_orientation == 5)
    {
        a = (drand(gen) - 0.5) * line_corr_angle * M_PI;
        k = tan(atan(k_first) + a);
        b = y1 - k * x1;
    }
}

/*----------------------------------------------------------------------------*/

void Linear_defects::add_first_line_correlated_case(int n,
                                                    complex<double> *x,
                                                    complex<double> *y,
                                                    double &k_first,
                                                    double &b_first) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    double x1 = 0.0;
    double y1 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;

    get_random_point_for_future_line(n, x, y, x1, y1);
    get_random_point_for_future_line(n, x, y, x2, y2);
    if(line_orientation == 3)
    {
        k_first = (y2 - y1) / (x2 - x1);
        b_first = y1 - k_first * x1;
    }
    if(line_orientation == 4)
    {
        y2 = y1 * (1.0 + 0.000000001 * (drand(gen) - 0.5));
        k_first = (y2 - y1) / (x2 - x1);
        b_first = y1 - k_first * x1;
    }
    if(line_orientation == 5)
    {
        x2 = x1 * (1.0 + 0.000000001 * (drand(gen) - 0.5));
        k_first = (y2 - y1) / (x2 - x1);
        b_first = y1 - k_first * x1;
    }
}






