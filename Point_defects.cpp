#include "Point_defects.h"

Point_defects::Point_defects() : io_data()
{
    get_input_point_defects_data(); 
}

/*----------------------------------------------------------------------------*/

Point_defects::~Point_defects()
{

}

/*----------------------------------------------------------------------------*/

void Point_defects::get_input_point_defects_data() 
{
    //gets all data form the point defects input file 
    string s;
    fstream input_file("./input/input_point_defects.txt");

    io_data.write_line_to_output_file("Input point defects data:\n");

    getline(input_file, s);
    potential = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Potential of LRDP  :" + to_string(potential));

    getline(input_file, s);
    d = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Range of potential :" + to_string(d));

    getline(input_file, s);
    rs_weak_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of  weak  RS  (1u) :" + to_string(rs_weak_conc) + " %");

    getline(input_file, s);
    rs_strong_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of strong RS (36u) :" + to_string(rs_strong_conc) + " %");

    getline(input_file, s);
    lrdp_weak_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of weak LRDP       :" + to_string(lrdp_weak_conc) + " %");

    getline(input_file, s);
    ci_positive_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of Coulomb CI(+)   :" + to_string(ci_positive_conc) + " %");

    getline(input_file, s);
    ci_negative_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of Coulomb CI(-)   :" + to_string(ci_negative_conc) + " %");

    getline(input_file, s);
    lrdp_strong_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of strong LRDP     :" + to_string(lrdp_strong_conc) + " %");

    getline(input_file, s);
    oxygen_graph_double_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of Oxygen, graphitic double subs.  :" + to_string(oxygen_graph_double_conc) + " %");

    getline(input_file, s);
    oxygen_graph_double2_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of Oxygen, graphitic double2 subs. :" + to_string(oxygen_graph_double2_conc) + " %");

    getline(input_file, s);
    oxygen_graph_double3_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of Oxygen, graphitic double3 subs. :" + to_string(oxygen_graph_double3_conc) + " %");

    getline(input_file, s);
    oxygen_monomeric_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of Oxygen, monomeric.              :" + to_string(oxygen_monomeric_conc) + " %");

    getline(input_file, s);
    oxygen_dimeric_conc = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file("Concentration of Oxygen, dimeric.                :" + to_string(oxygen_dimeric_conc) + " %");
    io_data.write_line_to_output_file("____________________\n"); 

    input_file.close();
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_point_defects(int n,
				      int n_x,
				      int n_y,
				      complex<double> *h_1,
			              complex<double> *h_2_u,
			              complex<double> *h_2_d,
			              complex<double> *h_3_u,
			              complex<double> *h_3_d,
				      complex<double> *x,
				      complex<double> *y) 
{
    if(rs_weak_conc     != 0.0) {add_rs_weak     (n, h_1);}
    if(rs_strong_conc   != 0.0) {add_rs_strong   (n, h_1);}
    if(lrdp_weak_conc   != 0.0) {add_lrdp_weak   (n, h_1, x, y);}
    if(ci_positive_conc != 0.0) {add_ci          ( 1.0, n, n_x, h_1, x, y);}
    if(ci_negative_conc != 0.0) {add_ci          (-1.0, n, n_x, h_1, x, y);}
    if(lrdp_strong_conc != 0.0) {add_lrdp_strong (n, h_1, x, y);}
    if(oxygen_graph_double_conc  != 0.0) {add_oxygen_graphitic_double  (n, n_x, n_y, h_1, h_2_u, h_2_d, h_3_u, h_3_d, x, y);}
    if(oxygen_graph_double2_conc != 0.0) {add_oxygen_graphitic_double2 (n, n_x, n_y, h_1, x, y);}
    if(oxygen_graph_double3_conc != 0.0) {add_oxygen_graphitic_double3 (n, n_x, n_y, h_1, x, y);}
    if(oxygen_monomeric_conc     != 0.0) {add_oxygen_monomeric         (n, n_x, n_y, h_1, h_2_u, h_2_d, h_3_u, h_3_d, x, y);}
    if(oxygen_dimeric_conc       != 0.0) {add_oxygen_dimeric           (n, n_x, n_y, h_1, h_2_u, h_2_d, h_3_u, h_3_d, x, y);}
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_rs_weak(int n,
				complex<double> *h_1) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    potential = 1.0;
    int rand_number = 0;

    for(int i = 0 ; i < int(n * rs_weak_conc / 100.0) ; i++)
    {
        rand_number = int(n * drand(gen));	
        h_1[rand_number] = potential;
    }
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_rs_strong(int n,
				  complex<double> *h_1) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    potential = 36.0;
    int rand_number = 0;

    for(int i = 0 ; i < int(n * rs_strong_conc / 100.0) ; i++)
    {
        rand_number = int(n * drand(gen));	
        h_1[rand_number] = potential;
    }
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_lrdp_weak(int n,
				  complex<double> *h_1,
				  complex<double> *x,
				  complex<double> *y) 
{
    int rand_number = 0;
    double delta_u = 1.0;
    double r = 0.0;
    potential = 1.0;
    d = 5.0;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    for(int i = 0 ; i < int(n * lrdp_weak_conc / 100.0) ; i++)
    {
        rand_number = int(n * drand(gen));	
        potential = 2.0 * delta_u * drand(gen) - delta_u;
        for(int j = 0 ; j < n ; j++)
        {
            r = get_ditance(real(x[rand_number]), 
                            real(x[j]), 
                            real(y[rand_number]), 
                            real(y[j]),
                            50.0);
            h_1[j] += potential * exp( - r * r / (2 * d * d));
        }
    }
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_ci(double ci_sign,
			   int n,
			   int n_x,
			   complex<double> *h_1,
			   complex<double> *x,
			   complex<double> *y) 
{
    /**
     * eps - dielectric permittivity of the SiO2
     * charged defects are located in a middle of hexagons
     * In x_ci and y_ci conserves information about position of ci centers
     *  
     * scaling_factor - all constants in eV
     *   
     * height_z - distance between grapehen and substrate
     */

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    double height_z   = 3.0;
    double eps        = 3.9;
    double eps_zero   = 8.85 * 1e-12;
    double e_electron = 1.6  * 1e-19;
    double a_lattice  = 1.42 * 1e-10;

    double scaling_factor = ci_sign / (4.0 * M_PI * eps * eps_zero / a_lattice);

    double amount_of_ci = 0.0;
    if(ci_sign > 0.0) {amount_of_ci = ci_positive_conc * n / 100.0;}
    if(ci_sign < 0.0) {amount_of_ci = ci_negative_conc * n / 100.0;}

    double x_ci = 0.0;
    double y_ci = 0.0;

    potential = 1.0;
    int rand_number = 0;

    double delta_u = 1.0;
    d = 5.0;

    double r = 0.0;

    for(int i = 0 ; i < amount_of_ci ; i++)
    {
        rand_number = int(n * drand(gen));	

        if((rand_number % 2 == 0) && (int(rand_number / n_x) % 2 == 0))
        {
            x_ci = real(x[rand_number]) + sqrt(3.0) / 2.0;
            y_ci = real(y[rand_number]) + 0.5;
        }
            if((rand_number % 2 != 0) && (int(rand_number / n_x) % 2 == 0))
        {
            x_ci = real(x[rand_number]);
            y_ci = real(y[rand_number]) + 1.0;
        }
        if((rand_number % 2 == 0) && (int(rand_number / n_x) % 2 != 0))
        {
            x_ci = real(x[rand_number]); 
            y_ci = real(y[rand_number]) + 1.0;
        }
        if((rand_number % 2 != 0) && (int(rand_number / n_x) % 2 != 0))
        {
            x_ci = real(x[rand_number]) + sqrt(3.0) / 2.0;
            y_ci = real(y[rand_number]) + 0.5;
        }				

        potential = 2.0 * delta_u * drand(gen) - delta_u;
        for(int j = 0 ; j < n ; j++)
        {
            r = get_ditance_for_ci(x_ci, 
                                   real(x[j]), 
                                   y_ci, 
                                   real(y[j]),
                                   height_z);
            h_1[j] += scaling_factor / r;
        }
    }
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_lrdp_strong(int n,
				    complex<double> *h_1,
				    complex<double> *x,
				    complex<double> *y) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    int rand_number = 0;
                  
    double r = 0.0;

    for(int i = 0 ; i < int(n * lrdp_strong_conc / 100.0) ; i++)
    {
        rand_number = int(n * drand(gen));	
        for(int j = 0 ; j < n ; j++)
        {
            r = get_ditance(real(x[rand_number]), 
                            real(x[j]), 
                            real(y[rand_number]), 
                            real(y[j]),
                            50.0);
            h_1[j] += potential * exp( - r * r / (2 * d * d));
        }
        //h_1[rand_number] = potential; 
    }
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_oxygen_graphitic_double(int n,
				   	        int n_x,
					        int n_y,
					        complex<double> *h_1,
					        complex<double> *h_2_u,
					        complex<double> *h_2_d,
					        complex<double> *h_3_u,
					        complex<double> *h_3_d,
					        complex<double> *x,
					        complex<double> *y) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);
    
    int n_n = int(n * oxygen_graph_double_conc / 100.0);

    int rand_number_1 = 0;
    int rand_number_2 = 0;

    double r = 0.0;
    int coin = 0;
    int top_down = 0; //top - 0, down - 1

    int *imps = new int[n_n];
    bzero(imps, sizeof(int) * n_n);

    for(int i = 0 ; i < int(n_n/2) ; i++)
    {
        rand_number_1 = get_random_site_not_close_to_borders(n, n_x, n_y);
        top_down = get_top_or_down_site(n, n_x, n_y, rand_number_1);
        coin = int(3 * drand(gen));

        if(coin == 0) 
        {
            rand_number_2 = rand_number_1 - 1;
            //_2_u[rand_number_1] = 0.0;
            //h_2_d[rand_number_1] = 0.0;
        }
        if(coin == 1) 
        {
            rand_number_2 = rand_number_1 + 1;
            //h_2_u[rand_number_1] = 0.0;
            //h_2_d[rand_number_1] = 0.0;
        }
        if ((coin == 2) && (top_down == 0)) 
        {
            rand_number_2 = rand_number_1 + n_x;
            //h_3_u[rand_number_1] = 0.0;
            //h_3_d[rand_number_1] = 0.0;
        }
        if ((coin == 2) && (top_down == 1)) 
        {
            rand_number_2 = rand_number_1 - n_x;
            //h_3_u[rand_number_1] = 0.0;
            //h_3_d[rand_number_1] = 0.0;
        }

        imps[i]         = rand_number_1;
        imps[i + n_n/2] = rand_number_2;
    }

    for(int i = 0 ; i < n_n ; i++)
    {
        for(int j = 0 ; j < n ; j++)
        {
            r = get_ditance(real(x[imps[i]]), 
                            real(x[j]), 
                            real(y[imps[i]]), 
                            real(y[j]),
                            50.0);
            if(r < 999998.0)
            {
                h_1[j] += potential * exp( - r * r / (2 * d * d));
            }
        }
    }

    if (imps != NULL) delete [] imps;    
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_oxygen_graphitic_double2(int n,
				   	         int n_x,
					         int n_y,
					         complex<double> *h_1,
					         complex<double> *x,
					         complex<double> *y) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);
    
    int n_n = int(n * oxygen_graph_double2_conc / 100.0);

    int rand_number_1 = 0;
    int rand_number_2 = 0;
    int rand_number_3 = 0;

    double r = 0.0;
    int coin = 0;
    int coin2 = 0;
    int top_down = 0; //top - 0, down - 1

    int *imps = new int[n_n];
    bzero(imps, sizeof(int) * n_n);

    for(int i = 0 ; i < int(n_n / 2) ; i++)
    {
        rand_number_1 = get_random_site_not_close_to_borders(n, n_x, n_y);
        top_down = get_top_or_down_site(n, n_x, n_y, rand_number_1);
        coin  = int(3 * drand(gen));
        coin2 = int(2 * drand(gen));

        if(coin == 0) 
        {
            rand_number_2 = rand_number_1 - 1;
            if(coin2 == 0) {rand_number_3 = rand_number_2 - 1;}    
            if(coin2 == 1) 
            {
                if(top_down == 0) {rand_number_3 = rand_number_2 + n_x;}
                if(top_down == 1) {rand_number_3 = rand_number_2 - n_x;}
            }                     
        }
        if(coin == 1)
        {
            rand_number_2 = rand_number_1 + 1;
            if(coin2 == 0) {rand_number_3 = rand_number_2 + 1;}    
            if(coin2 == 1) 
            {
                if(top_down == 0) {rand_number_3 = rand_number_2 + n_x;}
                if(top_down == 1) {rand_number_3 = rand_number_2 - n_x;}
            }                     
        }

        if ((coin == 2) && (top_down == 0)) 
        {
            rand_number_2 = rand_number_1 + n_x;
            if(coin2 == 0) {rand_number_3 = rand_number_2 - 1;}    
            if(coin2 == 1) {rand_number_3 = rand_number_2 + 1;}   
        }

        if ((coin == 2) && (top_down == 1)) 
        {
            rand_number_2 = rand_number_1 - n_x;
            if(coin2 == 0) {rand_number_3 = rand_number_2 - 1;}    
            if(coin2 == 1) {rand_number_3 = rand_number_2 + 1;}   
        }

        imps[i] = rand_number_1;
        imps[i + int(n_n / 2)] = rand_number_3;
    }

    for(int i = 0 ; i < n_n ; i++)
    {
        for(int j = 0 ; j < n ; j++)
        {
            r = get_ditance(real(x[imps[i]]), 
                            real(x[j]), 
                            real(y[imps[i]]), 
                            real(y[j]),
                            50.0);
            
            if(r < 999998.0)
            {
                h_1[j] += potential * exp( - r * r / (2 * d * d));
            }
            
        }
    }

    if (imps != NULL) delete [] imps;
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_oxygen_graphitic_double3(int n,
				   	         int n_x,
					         int n_y,
					         complex<double> *h_1,
					         complex<double> *x,
					         complex<double> *y) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);
    
    int n_n = int(n * oxygen_graph_double3_conc / 100.0);

    int rand_number_1 = 0;
    int rand_number_2 = 0;

    double r = 0.0;
    int coin = 0;
    int coin2 = 0;
    int top_down = 0; //top - 0, down - 1

    int *imps = new int[n_n];
    bzero(imps, sizeof(int) * n_n);

    for(int i = 0 ; i < int(n_n / 2) ; i++)
    {
        rand_number_1 = get_random_site_not_close_to_borders(n, n_x, n_y);
        top_down = get_top_or_down_site(n, n_x, n_y, rand_number_1);
        coin  = int(3 * drand(gen));
        coin2 = int(2 * drand(gen));

        if ((top_down == 0) && (coin == 0)) {rand_number_2 = rand_number_1 + n_x - 2;}
        if ((top_down == 1) && (coin == 0)) {rand_number_2 = rand_number_1 - n_x - 2;}
        if ((top_down == 0) && (coin == 1)) {rand_number_2 = rand_number_1 + n_x + 2;}
        if ((top_down == 1) && (coin == 1)) {rand_number_2 = rand_number_1 - n_x + 2;}
        if ((top_down == 0) && (coin == 2)) {rand_number_2 = rand_number_1 - n_x;}
        if ((top_down == 1) && (coin == 2)) {rand_number_2 = rand_number_1 + n_x;}

        imps[i] = rand_number_1;
        imps[i + int(n_n / 2)] = rand_number_2;
    }

    for(int i = 0 ; i < n_n ; i++)
    {
        for(int j = 0 ; j < n ; j++)
        {
            r = get_ditance(real(x[imps[i]]), 
                            real(x[j]), 
                            real(y[imps[i]]), 
                            real(y[j]),
                            50.0);
            
            if(r < 999998.0)
            {
                h_1[j] += potential * exp( - r * r / (2 * d * d));
            }
            
        }
        //h_1[imps[i]] = potential; 
    }

    if (imps != NULL) delete [] imps;
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_oxygen_monomeric(int n,
				   	 int n_x,
					 int n_y,
					 complex<double> *h_1,
					 complex<double> *h_2_u,
					 complex<double> *h_2_d,
					 complex<double> *h_3_u,
					 complex<double> *h_3_d,
					 complex<double> *x,
					 complex<double> *y) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);
    
    int n_n = int(n * oxygen_monomeric_conc / 100.0);

    int rand_number_1 = 0;
    int rand_number_2 = 0;

    double r = 0.0;
    int coin = 0;
    int top_down = 0; //top - 0, down - 1

    int *imps = new int[n_n];
    int *vacs = new int[n_n];
    bzero(imps, sizeof(int) * n_n);
    bzero(vacs, sizeof(int) * n_n);

    for(int i = 0 ; i < n_n ; i++)
    {
        rand_number_1 = get_random_site_not_close_to_borders(n, n_x, n_y);
        top_down = get_top_or_down_site(n, n_x, n_y, rand_number_1);
        coin = int(3 * drand(gen));

        if(coin == 0) {rand_number_2 = rand_number_1 - 1;}
        if(coin == 1) {rand_number_2 = rand_number_1 + 1;}

        if ((coin == 2) && (top_down == 0)) {rand_number_2 = rand_number_1 + n_x;}
        if ((coin == 2) && (top_down == 1)) {rand_number_2 = rand_number_1 - n_x;}

        imps[i]         = rand_number_1;
        vacs[i] = rand_number_2;
    }

    for(int i = 0 ; i < n_n ; i++)
    {
        for(int j = 0 ; j < n ; j++)
        {
            r = get_ditance(real(x[imps[i]]), 
                            real(x[j]), 
                            real(y[imps[i]]), 
                            real(y[j]),
                            50.0);
            
            if(r < 999998.0)
            {
                h_1[j] += potential * exp( - r * r / (2 * d * d));
            }
            
        }
    }

    for(int i = 0 ; i < n_n ; i++)
    {
        h_1  [vacs[i]] = 0.0;

        h_2_u[vacs[i]    ] = 0.0;
        h_2_u[vacs[i] - 1] = 0.0;
        h_2_d[vacs[i]    ] = 0.0;
        h_2_d[vacs[i] + 1] = 0.0;

        h_3_u[vacs[i]      ] = 0.0;
        h_3_u[vacs[i] - n_x] = 0.0;
        h_3_d[vacs[i]      ] = 0.0;
        h_3_d[vacs[i] + n_x] = 0.0;
    }

    if (imps != NULL) delete [] imps;
    if (vacs != NULL) delete [] vacs;   
}

/*----------------------------------------------------------------------------*/

void Point_defects::add_oxygen_dimeric(int n,
				   	 int n_x,
					 int n_y,
					 complex<double> *h_1,
					 complex<double> *h_2_u,
					 complex<double> *h_2_d,
					 complex<double> *h_3_u,
					 complex<double> *h_3_d,
					 complex<double> *x,
					 complex<double> *y) 
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);
    
    int n_n = int(n * oxygen_dimeric_conc / 100.0);

    int rand_number_1 = 0;
    int rand_number_2 = 0;
    int rand_number_3 = 0;

    double r = 0.0;
    int coin = 0;
    int coin2 = 0;
    int top_down = 0; //top - 0, down - 1

    int *imps = new int[n_n];
    int *vacs = new int[int(n_n / 2)];
    bzero(imps, sizeof(int) * n_n);
    bzero(vacs, sizeof(int) * int(n_n / 2));

    for(int i = 0 ; i < int(n_n / 2) ; i++)
    {
        rand_number_1 = get_random_site_not_close_to_borders(n, n_x, n_y);
        top_down = get_top_or_down_site(n, n_x, n_y, rand_number_1);
        coin  = int(3 * drand(gen));
        coin2 = int(2 * drand(gen));

        if(coin == 0) 
        {
            rand_number_2 = rand_number_1 - 1;
            if(coin2 == 0) {rand_number_3 = rand_number_2 - 1;}    
            if(coin2 == 1) 
            {
                if(top_down == 0) {rand_number_3 = rand_number_2 + n_x;}
                if(top_down == 1) {rand_number_3 = rand_number_2 - n_x;}
            }                     
        }
        if(coin == 1)
        {
            rand_number_2 = rand_number_1 + 1;
            if(coin2 == 0) {rand_number_3 = rand_number_2 + 1;}    
            if(coin2 == 1) 
            {
                if(top_down == 0) {rand_number_3 = rand_number_2 + n_x;}
                if(top_down == 1) {rand_number_3 = rand_number_2 - n_x;}
            }                     
        }


        if ((coin == 2) && (top_down == 0)) 
        {
            rand_number_2 = rand_number_1 + n_x;
            if(coin2 == 0) {rand_number_3 = rand_number_2 - 1;}    
            if(coin2 == 1) {rand_number_3 = rand_number_2 + 1;}   
        }

        if ((coin == 2) && (top_down == 1)) 
        {
            rand_number_2 = rand_number_1 - n_x;
            if(coin2 == 0) {rand_number_3 = rand_number_2 - 1;}    
            if(coin2 == 1) {rand_number_3 = rand_number_2 + 1;}   
        }

        imps[i] = rand_number_1;
        imps[i + int(n_n / 2)] = rand_number_3;
        vacs[i] = rand_number_2;
    }

    for(int i = 0 ; i < n_n ; i++)
    {
        for(int j = 0 ; j < n ; j++)
        {
            r = get_ditance(real(x[imps[i]]), 
                            real(x[j]), 
                            real(y[imps[i]]), 
                            real(y[j]),
                            50.0);
            
            if(r < 999998.0)
            {
                h_1[j] += potential * exp( - r * r / (2 * d * d));
            }
            
        }
        //h_1[imps[i]] = potential; 
    }

    for(int i = 0 ; i < int(n_n / 2) ; i++)
    {
        h_1  [vacs[i]      ] = 0.0;
        h_2_u[vacs[i]      ] = 0.0;
        h_2_u[vacs[i] - 1  ] = 0.0;
        h_2_d[vacs[i]      ] = 0.0;
        h_2_d[vacs[i] + 1  ] = 0.0;
        h_3_u[vacs[i]      ] = 0.0;
        h_3_u[vacs[i] - n_x] = 0.0;
        h_3_d[vacs[i]      ] = 0.0;
        h_3_d[vacs[i] + n_x] = 0.0;
    }

    if (imps != NULL) delete [] imps;
    if (vacs != NULL) delete [] vacs;
    
}

/*----------------------------------------------------------------------------*/

double Point_defects::get_ditance(double x1,
				  double x2,
				  double y1,
				  double y2,
				  double cutoff_r) 
{
    double distance = 0.0;
    double dx = x1 - x2;
    double dy = y1 - y2;

    if(dx < 0.0 ) {dx = -dx;}
    if(dx > cutoff_r) {return 999999.0;}
 
    if(dy < 0.0 ) {dy = -dy;}
    if(dy > cutoff_r) {return 999999.0;}

    distance = sqrt(dx * dx + dy * dy);

    return distance;
}

/*----------------------------------------------------------------------------*/

double Point_defects::get_ditance_for_ci(double x1,
				   	 double x2,
				   	 double y1,
				   	 double y2,
				   	 double height)
{
    double distance = 0.0;

    double dx = x1 - x2;
    double dy = y1 - y2;

    distance = sqrt(dx * dx + dy * dy + height * height);

    return distance;
}

/*----------------------------------------------------------------------------*/

double Point_defects::get_top_or_down_site(int n,
				   	   int n_x,
					   int n_y,
                                           int n_site)
{
    int td = 0; // 0 - up, 1 - down
    int n_y_site = int (n_site / n_x);
    int n_x_site = n_site - n_x * n_y_site;
    
    if (n_y_site % 2 == 0)
    {
        if (n_x_site % 2 == 0) {td = 0;}
        if (n_x_site % 2 != 0) {td = 1;}
    }
    if (n_y_site % 2 != 0)
    {
        if (n_x_site % 2 == 0) {td = 1;}
        if (n_x_site % 2 != 0) {td = 0;}
    }

    return td;
}

/*----------------------------------------------------------------------------*/

double Point_defects::get_random_site_not_close_to_borders(int n,
				   	                   int n_x,
					                   int n_y)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);


    int r_n  = 0;
    int r_nx = 0;
    int r_ny = 0;
    
    while(1)
    {
        r_n = int(n * drand(gen));
        r_ny = r_n / n_x;
        r_nx = r_n - r_ny * n_x;
        if ((n_y - r_ny > 5)
                && (r_ny       > 5)
                && (n_x - r_nx > 5)
                && (r_nx       > 5)) {break;}
    }

    return r_n;
}




