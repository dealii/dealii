// ======================================================
// MESH CONVERSION TOOL
//
// This program is distributed under the terms of the
// GNU GPL v2.0 license.
//
// Author: Jean-Paul Pelteret
//         jppelteret.uct@gmail.com
//         modified by: Timo Heister, heister@clemson.edu
// ======================================================

#include "MeshConversion.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

//--------------------------------------------------------------------------------

void display_help (void) {
	std::cout   << "The first input argument is the spatial dimension of the input file, "
		    << std::endl 
		    << "the second is the path to the read-in Abaqus .inp file, "
		    << std::endl 
		    << "and the third is the name of the file to which you wish to write the output AVS .ucd file to."
		    << std::endl 
		    << std::endl
		    << "Correct program usage is: " 
		    << std::endl
		    << "      './convert_mesh <spatial_dimension> /path/to/input_file.inp /path/to/output_file.ucd'"
		    << std::endl
		    << "e.g.  './convert_mesh 3 mesh/3d/test_in.inp mesh/3d/test_out.ucd'"
		    << std::endl
		    << "NOTE: New Abaqus files outputted by Cubit 12.0 and later MUST be generated with the 'NOT CUBIT ID's' option selected. "
		    << std::endl;
}

int main (int argc, char* argv[]) {
    if (argc != 4)
	{
	  display_help();
	  return 0;
	}
	
	try
	{
		const unsigned int dimension = atoi (argv[1]);
		std::cout << "Dimension: " << dimension << std::endl;
		std::string input_file = argv[2];
		std::string output_file = argv[3];
		
		MeshConversion mesh (dimension, input_file, output_file);
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
	
	return 0;
}

//--------------------------------------------------------------------------------
