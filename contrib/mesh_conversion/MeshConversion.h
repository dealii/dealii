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

#ifndef MESH_CONVERSION_H
#define MESH_CONVERSION_H

//--------------------------------------------------------------------------------

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

class MeshConversion
{
	public:
		MeshConversion (const unsigned int dimension, const std::string input_file, const std::string output_file);
		~MeshConversion ();
		bool convert_mesh ();
		
	private:
		void greeting ();

		bool read_in_abaqus ();
		bool write_out_avs_ucd ();
		
		std::vector <double> get_global_node_numbers(const int face_cell_no, const int face_cell_face_no);
		
		const double 				tolerance;
		const unsigned int 			dimension;
		
		unsigned int 				node_per_face;
		unsigned int 				node_per_cell;
		unsigned int 				face_per_cell;
		
		const std::string 			input_file_name;
		const std::string 			output_file_name;
		
		// NL: Stored as [ global node-id (int), x-coord, y-coord, z-coord ]
		std::vector< std::vector<double> > 	node_list; 
		// CL: Stored as [ material-id (int), node1, node2, node3, node4, node5, node6, node7, node8 ]	
		std::vector< std::vector<double> > 	cell_list; 
		// FL: Stored as [ sideset-id (int), node1, node2, node3, node4 ]	
		std::vector< std::vector<double> > 	face_list; 	
};

#endif // MESH_CONVERSION_H
