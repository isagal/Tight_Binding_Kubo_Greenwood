#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <sys/stat.h>


using namespace std;

#ifndef IO_data_H

class IO_data
{

public:
    IO_data();
    ~IO_data();
    
    void create_or_clean_output_directory();
    
    int get_int_value_from_input_line(string s);
    
    double get_double_value_from_input_line(string s);
    
    void write_line_to_output_file(string s);
    
    string get_string_value_from_input_line(string s);
};

#define IO_data_H
#endif // !IO_data_H




