#include "KuboGreenwood.h"
#include "IO_data.h"


int main()
{
    IO_data io_data;
    io_data.create_or_clean_output_directory();

    KuboGreenwood kb;
    kb.calculate_kb();

    return 0;
}
