#include "Psy.h"

Psy::Psy(Hamiltonian &h) : io_data()
{
    psy_null     = NULL;
    psy_null_dos = NULL;	
    psy_x_null   = NULL;

    psy          = NULL;
    psy_x        = NULL;
    x_psy        = NULL;
    x_psy_psy_x  = NULL;

    psy_null     = new complex<double>[h.n];
    psy_null_dos = new complex<double>[h.n];
    psy_x_null   = new complex<double>[h.n];

    psy          = new complex<double>[h.n];
    psy_x        = new complex<double>[h.n];
    x_psy        = new complex<double>[h.n];
    x_psy_psy_x  = new complex<double>[h.n];

    drop_all_arrays(h);	
}

/*----------------------------------------------------------------------------*/

Psy::~Psy()
{
    if (psy_null     != NULL) {delete [] psy_null;    }
    if (psy_null_dos != NULL) {delete [] psy_null_dos;}
    if (psy_x_null   != NULL) {delete [] psy_x_null;  }
    if (psy          != NULL) {delete [] psy;         }
    if (psy_x        != NULL) {delete [] psy_x;       }
    if (x_psy        != NULL) {delete [] x_psy;       }
    if (x_psy_psy_x  != NULL) {delete [] x_psy_psy_x; }
}

/*----------------------------------------------------------------------------*/

void Psy::drop_all_arrays(Hamiltonian &h)
{
    bzero(psy_null,     sizeof(complex<double>) * h.n);
    bzero(psy_null_dos, sizeof(complex<double>) * h.n);
    bzero(psy_x_null,   sizeof(complex<double>) * h.n);
    bzero(psy,          sizeof(complex<double>) * h.n);
    bzero(psy_x,        sizeof(complex<double>) * h.n);
    bzero(x_psy,        sizeof(complex<double>) * h.n);
    bzero(x_psy_psy_x,  sizeof(complex<double>) * h.n);
}

/*----------------------------------------------------------------------------*/

void Psy::renew_psy_null_arrays(Hamiltonian &h)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> drand(0.0, 1.0);

    double norm = 0.0;
    //Psy null for propagation   
    for(int i = (int)((h.n_y - h.dn_y) / 2.0) ; i < (int)((h.n_y + h.dn_y) / 2.0) ; i++)
    {
        for(int j = (i * h.n_x + int((h.n_x - h.dn_x) / 2.0)) ; j < i * h.n_x + int((h.n_x + h.dn_x) / 2.0) ; j++)
        {
            psy_null[j] = exp(complex<double>(0.0, 1.0) * 2.0 * M_PI * drand(gen)); 
        }
    }       
    norm = 1.0 / sqrt(h.dn_x * h.dn_y);
    cblas_zscal(h.n, &norm, psy_null, 1);

    //Psy null for DOS  
    for(int i = (int)((h.n_y - h.dn_y_dos) / 2.0) ; i < (int)((h.n_y + h.dn_y_dos) / 2.0) ; i++)
    {
        for(int j = (i * h.n_x + int((h.n_x - h.dn_x_dos) / 2.0)) ; j < i * h.n_x + int((h.n_x + h.dn_x_dos) / 2.0) ; j++)
        {
            psy_null_dos[j] = exp(complex<double>(0.0, 1.0) * 2.0 * M_PI * drand(gen)); 		
        }
    }       
    norm = sqrt(h.dn_x_dos * h.dn_y_dos);
    for(int i = 0 ; i < h.n ; i++) {psy_null_dos[i] /= norm;}

    for(int i = 0 ; i < h.n ; i++)
    {
        psy_x_null[i] = psy_null[i] * h.x[i];
    }
	
}

/*----------------------------------------------------------------------------*/

void Psy::print_psy_quadr(Hamiltonian &h)
{
/*	ofstream fres("./output/psy_quadr.dat", ios_base::out);

	for(int i = 0 ; i < h.n ; i++)
	{
		if((i + 1) % h.n_x == 0) {fres << psy_null_dos[i] << endl;}
		else {fres << psy_null_dos[i] << " ";}
	}	*/
}

/*----------------------------------------------------------------------------*/

complex<double> Psy::get_c_n(int n_c_n, 
                             double u,
                             double eigen_min, 
                             double eigen_max, 
                             double time)
{
    double h_bar = 6.582119281515e-16;    //eV    
    double a = (eigen_max - eigen_min) * u * time / (2.0 * h_bar);
    double b = (eigen_max + eigen_min) * u * time / (2.0 * h_bar); 
    complex<double> one_j(0.0,1.0); 
    complex<double> c_n(0.0,0.0);
    c_n = complex<double>(2.0,0.0) * exp(- b * one_j) * pow((-one_j), n_c_n) *complex<double>(jn(n_c_n, a),0.0);
    //cout << c_n << endl;
    return c_n;
}

/*----------------------------------------------------------------------------*/

void Psy::get_psy(Hamiltonian &h, 
                  double time,
                  complex<double> *psy_in,
                  complex<double> *psy_out)
{
    fstream fres;
    fres.open("./output/norma_outfile.dat", fstream::app|fstream::out|fstream::ate|fstream::binary);
	
    int number_of_psy_iterations = 100000;

    double norma_start = 0.0;
    double norma       = 0.0;
    double norma1      = 100000.0;
    double norma2      = 0.0;

    complex<double> complex_two(2.0,0.0);

    complex<double> c_n(0.0,0.0);

    complex<double> *temp = new complex<double>[h.n];

    complex<double> *phy_1 = new complex<double>[h.n];
    complex<double> *phy_2 = new complex<double>[h.n];	
    complex<double> *phy_3 = new complex<double>[h.n];	

    bzero(temp,    sizeof(complex<double>) * h.n);
    bzero(psy_out, sizeof(complex<double>) * h.n);
    bzero(phy_1,   sizeof(complex<double>) * h.n);
    bzero(phy_2,   sizeof(complex<double>) * h.n);
    bzero(phy_3,   sizeof(complex<double>) * h.n);

    norma_start = get_norm2_of_psy(h, psy_in);
    fres << "Norma start: " << norma_start << endl;

    //First iteration. c_n have been calculated using number of iteration = 0 

    memcpy(phy_1, psy_in, sizeof(complex<double>) * h.n);
    c_n = get_c_n(0, h.u, h.emin, h.emax, time);
    memcpy(temp, phy_1, sizeof(complex<double>) * h.n);	
    cblas_zscal(h.n, &c_n, temp, 1);
    memcpy(psy_out, temp, sizeof(complex<double>) * h.n);	

    //Second iteration. c_n have been calculated using number of iteration = 1

    h.mult_matrix(phy_1, phy_2, h.diags_norm);
    c_n = get_c_n(1, h.u, h.emin, h.emax, time);
    memcpy(temp, phy_2, sizeof(complex<double>) * h.n);	
    cblas_zscal(h.n, &c_n, temp, 1);
    h.add_array_x_to_y (psy_out, temp);

    //Enth iteration. c_n have been calculated using number_of_psy_iterations 
	
    for(int i = 2 ; i < number_of_psy_iterations ; i++)
    {
    	h.mult_matrix(phy_2, phy_3, h.diags_norm);
        cblas_zscal(h.n, &complex_two, phy_3, 1);
        h.subtract_array_x_form_y (phy_3, phy_1);
        c_n = get_c_n(i, h.u, h.emin, h.emax, time);	
    	memcpy(temp, phy_3, sizeof(complex<double>) * h.n);	
        cblas_zscal(h.n, &c_n, temp, 1);		
        h.add_array_x_to_y (psy_out, temp);		
    	
        memcpy(phy_1, phy_2, sizeof(complex<double>) * h.n);
    	memcpy(phy_2, phy_3, sizeof(complex<double>) * h.n);		

        norma = get_norm2_of_psy(h, psy_out);
		
        if (int(i % 50) == 0)
        {
            norma1 = norma2;
            norma2 = norma;
            fres << "    " << i + 1 << "    " << norma << endl;
        }
        if((fabs(norma2 - norma1) < 1e-15 * norma_start)) 
        {
            if (temp  != NULL) delete [] temp;
            if (phy_1 != NULL) delete [] phy_1;
            if (phy_2 != NULL) delete [] phy_2;
            if (phy_3 != NULL) delete [] phy_3;
            return;
        }     
    }
    io_data.write_line_to_output_file("Norma did not converge correctly"); 	

    if (temp  != NULL) delete [] temp;
    if (phy_1 != NULL) delete [] phy_1;
    if (phy_2 != NULL) delete [] phy_2;
    if (phy_3 != NULL) delete [] phy_3;
}

/*----------------------------------------------------------------------------*/

double Psy::get_norm2_of_psy(Hamiltonian &h,
                             complex<double> *v)
{
    complex<double> dnorm(0.0,0.0);
    int incx = 1, incy = 1;
    zdotc(&dnorm, &h.n, v, &incx, v, &incy);
    return dnorm.real();
}

/*----------------------------------------------------------------------------*/

void Psy::get_x_psy_psy_x(Hamiltonian &h)
{
    for(int i = 0 ; i < h.n ; i++)
    {
        x_psy[i] = psy[i] * h.x[i];
    }	
    memcpy(x_psy_psy_x, x_psy, sizeof(complex<double>) * h.n); 		
    h.subtract_array_x_form_y (x_psy_psy_x, psy_x);
}



/*
double Psy::normalize()
{
	double dnorm = norm();
	kbVec c(LongiN, TransN);
	int incx = 1, incy = 1;
	complex<double> a(1.0/sqrt(dnorm),0.0);
	zaxpy(&size, &a, v, &incx, c.v, &incy);
	zcopy(&size, c.v, &incx, v, &incy);
	return dnorm;
}*/

