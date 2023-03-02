#include "IO_data.h"

IO_data::IO_data() 
{

}

IO_data::~IO_data()
{

}

void IO_data::create_or_clean_output_directory() 
{
    //creates or cleans output directory
    struct stat st;
    if(stat("./output",&st) == 0)
    {
        remove("./output/output.txt");
        remove("./output/h_1.dat");
        remove("./output/psy_quadr.dat");
        remove("./output/dos.dat");
        remove("./output/continued_fraction_test.dat");
        remove("./output/norma_outfile.dat");
        remove("./output/diffusivity.dat");
        remove("./output/conductivity.dat");
        remove("./output/conductivity_in_simenses.dat");
    }       
    else
    {
        mkdir("output", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

}

int IO_data::get_int_value_from_input_line(string s) 
{
    //Returns int value at the end of the line s with a delimeter ": "
    return stoi(get_string_value_from_input_line(s));
}

double IO_data::get_double_value_from_input_line(string s) 
{
    //Returns double value at the end of the line s with a delimeter ": "
    return stod(get_string_value_from_input_line(s));
}

void IO_data::write_line_to_output_file(string s) 
{
    //writes string s to the output.txt
    fstream output_file;
    output_file.open("./output/output.txt", fstream::app|fstream::out|fstream::ate|fstream::binary);
    output_file << s << endl;
    output_file.close();
}

string IO_data::get_string_value_from_input_line(string s)
{
    //returns value from the input line accordingly to a delimiter ": "
    string delimiter = ": ";
    s.erase(0, s.find(delimiter) + delimiter.length());
    //cout << "_" << s << "_" << endl;
    return s;
}

