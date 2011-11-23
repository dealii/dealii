#include <fstream>
#include <iostream>

int main ()
{
    std::ofstream logfile("double_bug/output");

    logfile <<  -0.00976562 << std::endl;

    return 0;
}
