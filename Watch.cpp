#include "Watch.h"

Watch::Watch() 
{
    srand (time(NULL));
    time_spent_begin = clock();
    time_spent_end   = clock();
}

Watch::~Watch()
{

}

void Watch::get_and_print_work_time(IO_data io_data) 
{
    io_data.write_line_to_output_file("Total work time: " + get_string_work_time() );  
}

string Watch::get_string_work_time() 
{
    string s(""); 
    time_spent_end = clock();

    double work_time = (time_spent_end - time_spent_begin) / CLOCKS_PER_SEC;
    int hours   = work_time / 3600;
    int minutes = (work_time - hours * 3600) / 60;
    int seconds = work_time - hours * 3600 - minutes * 60;

    s = to_string(hours)   + " hours " 
      + to_string(minutes) + " min. " 
      + to_string(seconds) + " sec. ";

    return s;  
}


