#include "MeshConversion.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

//--------------------------------------------------------------------------------

void display_help (void) {
	std::cout   << "The first input argument is the spatial dimension of the input file, "
		    << std::endl 
		    << "the second is the type of input file (0 for Abaqus OLD, 1 for Abaqus NEW), "
		    << std::endl 
		    << "the third is the path to the read-in Abaqus .inp file, "
		    << std::endl 
		    << "and the fourth is the name of the file to which you wish to write the output AVS .ucd file to." 
		    << std::endl 
		    << std::endl
		    << "Correct program usage is: " 
		    << std::endl
		    << "      './convert_mesh <spatial_dimension> <input_file_type> /path/to/input_file.inp /path/to/output_file.ucd'" 
		    << std::endl
		    << "e.g.  './convert_mesh 3 0 mesh/3d/test_in.inp mesh/3d/test_out.ucd'" 
		    << std::endl
		    << "NOTE: Old Abaqus files are outputted by Cubit 11.1 and earlier. "
		    << std::endl
		    << "      New Abaqus files are outputted by Cubit 12.0 and later. They MUST be generated with the 'NOT CUBIT ID's' option selected. "
		    << std::endl;
}

int main (int argc, char* argv[]) {
    if (argc == 5)
    {
	try
	{
		const unsigned int dimension = atoi (argv[1]); //*argv[1] - '0';
		std::cout << "Dimension: " << dimension << std::endl;
		const unsigned int input_type = atoi (argv[2]); //*argv[2] - '0';
		std::string input_file = argv[3];
		std::string output_file = argv[4];
		
		MeshConversion mesh (dimension, input_type, input_file, output_file);
		const bool success = mesh.convert_mesh ();
		
		if (success == false) {
			std::cout << "Something has gone wrong with the conversion..." << std::endl;
			display_help();
		}
	}

	catch (...)
	{
	    std::cerr << std::endl
		      << "==============" << std::endl;
	    std::cerr << "UNKNOWN ERROR!" << std::endl
		      << "==============" << std::endl;

	    return 1;
	}
    }
    else if (argc == 2) {
    	if (argv[1] == "-h" || argv[1] == "-help")
    		display_help();
	else
    		display_help();
    }
    else {
	std::cerr   << "WRONG NUMBER OF INPUT ARGUMENTS:" << std::endl;
	display_help();		    
    }
	
	return 0;
}

//--------------------------------------------------------------------------------
