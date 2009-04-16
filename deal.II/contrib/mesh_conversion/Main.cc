#include "MeshConversion.h"
#include <iostream>

using namespace std;

//--------------------------------------------------------------------------------

int main (int argc, char* argv[]) {

    if (argc == 4)
    {
	try
	{
		const unsigned int dimension = *argv[1] - '0';
		std::cout << "Dimension: " << dimension << std::endl;
		std::string input_file = argv[2];
		std::string output_file = argv[3];
		
		MeshConversion mesh (dimension);

		const bool readin_successul = mesh.read_in_abaqus_inp(input_file);

		if (readin_successul == true)
			mesh.write_out_avs_ucd(output_file);
		else
			std::cerr << "ERROR: Input file not found." << std::endl;
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
    else {
	std::cerr   << "WRONG NUMBER OF INPUT ARGUMENTS:" << std::endl
		    << "The first input argument is the spatial dimension of the input file, the second is the path to the read-in Abaqus .inp file, and the third is the name of the file to which you wish to write the output AVS .ucd file to." << std::endl << std::endl
		    << "Correct program usage is: " << std::endl
		    << "      './convert_mesh <spatial_dimension> /path/to/input_file.inp /path/to/output_file.ucd'" << std::endl
		    << "e.g.  './convert_mesh 3 mesh/3d/test_in.inp mesh/3d/test_out.ucd'" << std::endl;
    }
	
	return 0;
}

//--------------------------------------------------------------------------------
