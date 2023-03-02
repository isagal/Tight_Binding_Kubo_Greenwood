#pragma once

#include <time.h>
#include <fstream>
#include <string>

#include "IO_data.h"


using namespace std;

#ifndef Watch_H

class Watch
{

public:
    Watch();
    ~Watch();

    IO_data io_data;

    double time_spent_begin, time_spent_end; 

    void get_and_print_work_time(IO_data io_data);
    string get_string_work_time(); 
};

#define Watch_H
#endif // !Watch_H




