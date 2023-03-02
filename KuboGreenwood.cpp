#include "KuboGreenwood.h"

KuboGreenwood::KuboGreenwood() : io_data(), ham() , psy(ham), watch()
{
    dos     = NULL;
    el_conc = NULL;

    get_input_kb_data();

    energy_start = -6.0;
    energy_end   =  6.0;

    dos       = new double[ham.n_energy_points];
    el_conc   = new double[ham.n_energy_points]; 
    diff      = new double[ham.n_energy_points];  
    cond      = new double[ham.n_energy_points]; 
    cond_in_s = new double[ham.n_energy_points]; 
}

/*----------------------------------------------------------------------------*/

KuboGreenwood::~KuboGreenwood()
{
    if (dos       != NULL) {delete [] dos;      }
    if (el_conc   != NULL) {delete [] el_conc;  }
    if (diff      != NULL) {delete [] diff;     }
    if (cond      != NULL) {delete [] cond;     }
    if (cond_in_s != NULL) {delete [] cond_in_s;}
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::drop_all_arrays()
{
    bzero(dos,       sizeof(double) * ham.n_energy_points);
    bzero(el_conc,   sizeof(double) * ham.n_energy_points);
    bzero(diff,      sizeof(double) * ham.n_energy_points);
    bzero(cond,      sizeof(double) * ham.n_energy_points);
    bzero(cond_in_s, sizeof(double) * ham.n_energy_points);
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::get_input_kb_data() 
{
    //gets all data form the Kubo-Greenwood input file 
    string s;
    fstream input_file("./input/input_kb.txt");

    io_data.write_line_to_output_file("Input KB data:\n");

    getline(input_file, s);
    flag_dos = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Calculation of DOS          (1 - yes, 0 - no) : ") + to_string(flag_dos));

    getline(input_file, s);
    flag_diff = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Calculation of diffusivity  (1 - yes, 0 - no) : ") + to_string(flag_diff));

    getline(input_file, s);
    flag_cond = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Calculation of conductivity (1 - yes, 0 - no) : ") + to_string(flag_cond));

    getline(input_file, s);
    diff_start_time = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Diffusivity starting time (in fs) : ") + to_string(diff_start_time));

    getline(input_file, s);
    diff_step_time = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Diffusivity time step     (in fs) : ") + to_string(diff_step_time));

    getline(input_file, s);
    diff_n_steps = io_data.get_int_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Diffusivity number of steps       : ") + to_string(diff_n_steps));

    getline(input_file, s);
    cond_time = io_data.get_double_value_from_input_line(s);
    io_data.write_line_to_output_file(string("Conductivity time (in fs) : ") + to_string(cond_time));

    io_data.write_line_to_output_file("____________________\n"); 

    diff_start_time /= 1e15;
    diff_step_time  /= 1e15;
    cond_time       /= 1e15;

    input_file.close();
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::calculate_kb()
{
    if(flag_dos  == 1) {calculate_dos()         ;}
    if(flag_diff == 1) {calculate_diffusivity() ;}
    if(flag_cond == 1) {calculate_conductivity();}
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::calculate_dos()
{
    string tot_time;

    drop_all_arrays();
    write_cont_frac_test();

    for(int i = 0 ; i < ham.n_realizations ; i++)
    {
        ham.drop_all_arrays();
        ham.get_lattice();
        ham.get_energy_array(energy_start, energy_end);
        psy.renew_psy_null_arrays(ham);
        ham.get_alpha_beta(psy.psy_null_dos);

        for(int j = 0 ; j < ham.n_energy_points ; j++)
        {
            dos[j] += calculate_dos_for_selected_energy(ham.energy_array[j]);
        }

        tot_time = watch.get_string_work_time(); 
        io_data.write_line_to_output_file(to_string(i) + "th realization has finished: " + tot_time);
    }
	
    for(int j = 0 ; j < ham.n_energy_points ; j++) {dos[j] /= ham.n_realizations;}

    print_dos();
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::calculate_diffusivity()
{
    double time = diff_start_time;
    double norma = 0.0;
    string tot_time;

    drop_all_arrays();
    write_cont_frac_test();

    for(int t = 0 ; t < diff_n_steps ; t++)
    {
        time = diff_start_time + t * diff_step_time;

        for(int i = 0 ; i < ham.n_realizations ; i++)
        {
            ham.drop_all_arrays();
            ham.get_lattice();
            ham.get_energy_array(energy_start, energy_end);
            psy.drop_all_arrays(ham);
            psy.renew_psy_null_arrays(ham);
		
            norma = 0.0;
            get_psy_propagation(time, norma);

            ham.get_alpha_beta(psy.psy_null_dos);
            for(int j = 0 ; j < ham.n_energy_points ; j++)
            {
                dos[j] += calculate_dos_for_selected_energy(ham.energy_array[j]);
            }
            ham.get_alpha_beta(psy.x_psy_psy_x);
            for(int j = 0 ; j < ham.n_energy_points ; j++)
            {
                diff[j] += calculate_diff_for_selected_energy_and_time(ham.energy_array[j], time, norma);
            }

            tot_time = watch.get_string_work_time(); 
            io_data.write_line_to_output_file(to_string(i) 
                    + "th realization for time " 
                    + to_string(time*1e15)
                    + " fs has finished: " 
                    + tot_time);
        }
        for(int j = 0 ; j < ham.n_energy_points ; j++)
        {
            diff[j] /= dos[j];
        }

        print_diff(t, time);
        drop_all_arrays();
    }
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::calculate_conductivity()
{
    double time = cond_time;
    double norma = 0.0;
    string tot_time;

    drop_all_arrays();
    write_cont_frac_test();

    for(int i = 0 ; i < ham.n_realizations ; i++)
    {
        ham.drop_all_arrays();
        ham.get_lattice();
        ham.get_energy_array(energy_start, energy_end);
        psy.drop_all_arrays(ham);
        psy.renew_psy_null_arrays(ham);
		
        norma = 0.0;
        get_psy_propagation(time, norma);

        ham.get_alpha_beta(psy.psy_null_dos);
        for(int j = 0 ; j < ham.n_energy_points ; j++)
        {
            dos[j] += calculate_dos_for_selected_energy(ham.energy_array[j]);
        }
        ham.get_alpha_beta(psy.x_psy_psy_x);
        for(int j = 0 ; j < ham.n_energy_points ; j++)
        {
            cond[j] += calculate_cond_for_selected_energy_and_time(ham.energy_array[j], time, norma);
        }

        tot_time = watch.get_string_work_time(); 
        io_data.write_line_to_output_file(to_string(i) 
                + "th realization for time " 
                + to_string(time*1e15)
                + " fs has finished: " 
                + tot_time);
    }

    for(int j = 0 ; j < ham.n_energy_points ; j++)
    {
        dos[j]  /= ham.n_realizations;
        cond[j] /= ham.n_realizations; 
        cond_in_s[j] = cond[j] * 7.748 / 100000.0 / 2.0;
    }
	
    calculate_el_conc();
    get_min_cond_in_s();
    print_cond();
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::get_min_cond_in_s()
{
    /**
     *Minimal conductivity in S is obtained within the energy range from -2.5 eV to 2.5 eV
     */

    double minimal_cond = 1000000.0;
    double minimal_cond_energy = -20.0;

    for(int i = 0 ; i < ham.n_energy_points ; i++)
    {
        if((ham.energy_array[i]* ham.u > -2.5) && (ham.energy_array[i]* ham.u < 2.5))
        {
            if(minimal_cond > cond[i])
            {
                minimal_cond_energy = ham.energy_array[i]; 
                minimal_cond = cond[i];
            }
        }
    }
                
    minimal_cond_energy *= ham.u;
    minimal_cond *= 7.748 / 100000.0 / 2.0;    
                                   
    io_data.write_line_to_output_file("\n\nMinimal conductivity in siemenses: " 
            + to_string(minimal_cond)
            + "S for energy " 
            + to_string(minimal_cond_energy)
            + "eV");
}
                          
/*----------------------------------------------------------------------------*/

double KuboGreenwood::calculate_dos_for_selected_energy(double energy)
{
    double temp_dos = 0.0;
    temp_dos = - 2.0 / M_PI * ham.get_continued_fraction(ham.cont_frac_n, energy);
    return temp_dos;
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::get_psy_propagation(double time, double &norma)
{
    complex<double> znorma(0.0,0.0);

    psy.get_psy(ham, time, psy.psy_null  , psy.psy  );
    psy.get_psy(ham, time, psy.psy_x_null, psy.psy_x); 
    psy.get_x_psy_psy_x(ham);

    norma = sqrt(psy.get_norm2_of_psy(ham, psy.x_psy_psy_x));
	
    znorma = complex<double> (norma,0.0);
    znorma = complex<double> (1.0,0.0) / znorma;
    cblas_zscal(ham.n, &znorma, psy.x_psy_psy_x, 1);
}

/*----------------------------------------------------------------------------*/

double KuboGreenwood::calculate_diff_for_selected_energy_and_time(double energy,
                                                                  double time,
                                                                  double &norma)
{
    double coef_radchenko = 0.75;
    double coef_spatial   = 2.0;	
    double temp_diff = 0.0;

    temp_diff  = -1.0 * norma * norma / (M_PI * time);
    temp_diff *= ham.get_continued_fraction(ham.cont_frac_n, energy);
    temp_diff *= 1e-15 * coef_radchenko / coef_spatial;

    return temp_diff;
}

/*----------------------------------------------------------------------------*/

double KuboGreenwood::calculate_cond_for_selected_energy_and_time(double energy,
                                                                  double time,
                                                                  double &norma)
{
    double coef_radchenko = 0.75;
    double coef_spatial   = 2.0;	
    double temp_cond = 0.0;

    temp_cond  = -1 * norma * norma * 6.6260695759e-34; 
    temp_cond /= (M_PI * time * 1.6 * 1e-19 * 3.0 * sqrt(3.0) / 2.0);
    temp_cond *= ham.get_continued_fraction(ham.cont_frac_n, energy);
    temp_cond *= coef_radchenko / coef_spatial;

    return temp_cond;
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::calculate_el_conc()
{
    for(int i = 0 ; i < (ham.n_energy_points - 1) ; i++)
    {
        el_conc[i + 1] = el_conc[i] + (dos[i] + (dos[i + 1] - dos[i])/2.0) * (energy_end - energy_start) / (ham.n_energy_points);

    }
    for(int i = 0 ; i < ham.n_energy_points ; i++)
    {
        el_conc[i] -= 1.0;
    }
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::print_dos()
{
    ofstream fres("./output/dos.dat", ios_base::out);
    fres << "energy dos" << endl;

    for(int i = 0 ; i < ham.n_energy_points ; i++)
    {
        fres << ham.energy_array[i] << " " << dos[i] << endl;
    }	
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::print_cond()
{
    ofstream fres1("./output/conductivity.dat", ios_base::out);
    fres1 << "energy dos el_conc cond" << endl;

    for(int i = 0 ; i < ham.n_energy_points ; i++)
    {
        fres1 << ham.energy_array[i] << " " 
                << dos[i]     << " " 
                << el_conc[i] << " " 
                << cond[i]    << endl;
    }	

    ofstream fres2("./output/conductivity_in_simenses.dat", ios_base::out);
    fres2 << "energy(eV) cond(S)" << endl;

    for(int i = 0 ; i < ham.n_energy_points ; i++)
    {
        fres2 << ham.energy_array[i] * ham.u << " "
                << el_conc[i] << " " 
                << el_conc[i] * 3.9e15 << " "  
                << cond_in_s[i] << endl;
    }
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::print_diff(int n_time_step, 
							   double time)
{
    fstream fres;
    fres.open("./output/diffusivity.dat", fstream::app|fstream::out|fstream::ate|fstream::binary);

    if(n_time_step == 0)
    {
        fres << "time(fs)" << " ";
        fres << ham.energy_array[(ham.n_energy_points-1)/2 - (ham.n_energy_points-1)/10] << " ";
        fres << ham.energy_array[(ham.n_energy_points-1)/2 - (ham.n_energy_points-1)/20] << " ";
        fres << ham.energy_array[(ham.n_energy_points-1)/2] << " ";
        fres << ham.energy_array[(ham.n_energy_points-1)/2 + (ham.n_energy_points-1)/20] << " ";
        fres << ham.energy_array[(ham.n_energy_points-1)/2 + (ham.n_energy_points-1)/10] << endl;
    }
    fres << time * 1e15 << " ";
    fres << diff[(ham.n_energy_points-1)/2 - (ham.n_energy_points-1)/10] << " ";
    fres << diff[(ham.n_energy_points-1)/2 - (ham.n_energy_points-1)/20] << " ";
    fres << diff[(ham.n_energy_points-1)/2] << " ";
    fres << diff[(ham.n_energy_points-1)/2 + (ham.n_energy_points-1)/20] << " ";
    fres << diff[(ham.n_energy_points-1)/2 + (ham.n_energy_points-1)/10] << endl;
}

/*----------------------------------------------------------------------------*/

void KuboGreenwood::write_cont_frac_test()
{
    /**
     *Method calculates and prints cont_frac values for convergence checking.
     */
    string tot_time;
    ofstream fres("./output/continued_fraction_test.dat", ios_base::out);

    int cfn = 0;
    int n_cf_en_points = 5;
    double temp_dos = 0.0;
    double *energy_array_con_frac = new double[n_cf_en_points];

    for(int j = 0 ; j < n_cf_en_points ; j++)
    {
        energy_array_con_frac[j] = -2 + 4 * j / (n_cf_en_points - 1);
    }

    ham.drop_all_arrays();
    ham.get_lattice();
    psy.renew_psy_null_arrays(ham);
    ham.get_alpha_beta(psy.psy_null_dos);

    fres << "n ";
    for(int i = 0 ; i < n_cf_en_points; i++) {fres << energy_array_con_frac[i] << " " ;}
    fres << endl;

    for(int i = 0 ; i < 11 ; i++)
    {
        cfn = (int)(2.0 + ham.cont_frac_n * i / 11.0);
        fres << cfn << " ";
        for(int j = 0 ; j < n_cf_en_points ; j++)
        {
            temp_dos = -2.0 / M_PI * ham.get_continued_fraction(cfn, energy_array_con_frac[j]);
            fres << temp_dos << " ";
        }
        fres << endl;
    }

    tot_time = watch.get_string_work_time(); 
    io_data.write_line_to_output_file("Continued fraction test have been completed: " + tot_time);

    if (energy_array_con_frac != NULL) delete [] energy_array_con_frac;
}






