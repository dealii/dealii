#include "MeshConversion.h"

//--------------------------------------------------------------------------------

MeshConversion::MeshConversion (const unsigned int dimension):
dimension (dimension),
tolerance (5e-16) // Used to offset Cubit tolerance error when outputting value close to zero
{
	greeting ();
	
	if (dimension == 2)
	{
		node_per_face = 2; 		// Line edge
		node_per_cell = 4;		// Quad face
		face_per_cell = 2*dimension;	// Lines per quad
	}
	else if (dimension == 3)
	{
		node_per_face = 4; 		// Quad face
		node_per_cell = 8; 		// Hex cell
		face_per_cell = 2*dimension; 	// Quads per hex
	}
	else
	{
		std::cout << "ERROR: Chosen spatial dimension is invalid!" << std::endl;
	}
}

//--------------------------------------------------------------------------------

MeshConversion::~MeshConversion ()
{

}

//--------------------------------------------------------------------------------

// ========================================================
// http://www.codeguru.com/forum/showthread.php?t=231054
// ========================================================
template <class T> bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

//--------------------------------------------------------------------------------

void MeshConversion::greeting () {
	std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;
	std::cout << "************************************************************************************" << std::endl;
	std::cout << "***                         FEM Mesh conversion tool                             ***" << std::endl;
	std::cout << "************************************************************************************" << std::endl;
	std::cout << "*** Author: Jean-Paul Pelteret                                                   ***" << std::endl;
	std::cout << "*** Date:   December 2008                                                        ***" << std::endl;
	std::cout << "*** References:                                                                  ***" << std::endl;
	std::cout << "***    http://www.dealii.org/developer/doxygen/deal.II/structGeometryInfo.html   ***" << std::endl;
	std::cout << "***    http://people.scs.fsu.edu/~burkardt/html/ucd_format.html                  ***" << std::endl;
	std::cout << "***    http://people.scs.fsu.edu/~burkardt/data/ucd/ucd.html                     ***" << std::endl;
	std::cout << "***    http://www.cprogramming.com/tutorial/string.html                          ***" << std::endl;
	std::cout << "***    http://www.codeguru.com/forum/showthread.php?t=231054                     ***" << std::endl;
	std::cout << "************************************************************************************" << std::endl;
	std::cout << "*** FEATURES:                                                                    ***" << std::endl;
	std::cout << "*** Read-in file types:                                                          ***" << std::endl;
	std::cout << "***      - Abaqus inp                                                            ***" << std::endl;
	std::cout << "*** Write-out file types:                                                        ***" << std::endl;
	std::cout << "***      - AVS UCD                                                               ***" << std::endl;
	std::cout << "************************************************************************************" << std::endl;
	std::cout << std::endl;
}

//--------------------------------------------------------------------------------

bool MeshConversion::read_in_abaqus_inp (const std::string filename) 
{
	input_file_name = filename;
	bool read_in_successful = true;
	std::cout << "Reading in ABAQUS .inp FILE: " << input_file_name << std::endl;
	
	std::ifstream input_stream;
	
	input_stream.open(filename.c_str(), std::ifstream::in);
	
	// Check to see if file exists
	if (input_stream.good() == true)
	{
		std::string buffer;
		std::vector < std::string > temp_data;
		
		// Loop over all the contents of the fibre data file and read it into buffer
		while (input_stream >> buffer)
		{
			for(int i = 0; i < buffer.length(); ++i)
				if (buffer[i] == ',')	// Get rid of the .inp file's useless comma's
				{
					buffer.erase(i);
					--i;
				}
				
			temp_data.push_back(buffer);
		}
		
		const bool print_input_deck = false;
		if (print_input_deck == true)
		{
			// =========== TEMP ===========
			std::string filename_out ("checks/input_stream.txt");
			std::ofstream output;
			output.open(filename_out.c_str());
			
			for (int ii = 0; ii < temp_data.size(); ii++ )
				output << temp_data[ii] << std::endl;
			
			output.close();
			// =========== TEMP ===========
		}
		
		for (int k = 0; k < temp_data.size(); ++k)
		{
			// ================================ NODES ===================================
			// ABAQUS formatting
			// 	*NODE
			// 	<global node number>, <x coordinate>, <y coordinate>, <z coordinate>
			// ==========================================================================

			if (temp_data[k] == "*NODE")
			{
				int j = 0;
				float temp_float;
				
				while (from_string<float> (temp_float, temp_data[k + 1 + j*(dimension+1+(dimension==2 ? 1 : 0))], std::dec) == true)
				{
					// Initilise storage variables
					std::vector <double> node (dimension+1);
					
					// Convert from string to the variable types
					for (int i = 0; i < dimension+1; ++i)
						from_string<double> (node[i], temp_data[k + 1 + j*(dimension+1+(dimension==2 ? 1 : 0)) + i], std::dec);
						
					// Add to the global node number vector
					node_list.push_back(node);
					
					++j;
				}
			}
			
			// ===================== 2D QUADS / 3D HEX CELLS =======================
			// ABAQUS formatting
			// 	*ELEMENT, TYPE=C3D8R, ELSET=EB< block-set number (material-id) >
			// 	<cell number>, <8 x node number>
			// =====================================================================

			else if (temp_data[k] == "*ELEMENT")
			{
				int j = 0;
				float temp_float;
				const int data_per_cell = node_per_cell + 1;
				
				// Get material id
				std::string material_id_line = temp_data[k+ 2];
				std::string material_id_temp;
				for (int ll = 8 /*Characters in "ELSET=EB" */; ll < material_id_line.length(); ++ll)
					material_id_temp += material_id_line[ll];
				
				int material_id = 0;
				for (int ll = 0; ll < material_id_temp.length(); ++ll)
					material_id += (material_id_temp[material_id_temp.length() - ll - 1] - 48 /* ASCII TRICKS */) * pow(10.0,ll);
				
				while (from_string<float> (temp_float, temp_data[k + 3 + j*(data_per_cell)], std::dec) == true)
				{
					// Initilise storage variables
					std::vector <double> cell (data_per_cell);
					
					// Material id
					cell[0] = material_id;
					
					// Convert from string to the variable types
					for (int i = 1; i < data_per_cell; ++i)
						from_string<double> (cell[i], temp_data[k + 3 + j*(data_per_cell) + i], std::dec);
						
					// Add to the global node number vector
					cell_list.push_back(cell);
					
					++j;
				}
			}
			
			// =========== 2D LINES / 3D BC_QUADS ===========
			// ABAQUS formatting
			// 	SURFACE, NAME=SS<side-set number>
			// 	<cell number>, <cell side number>
			// ==============================================
			
			else if (temp_data[k] == "*SURFACE" )
			{
				float temp_float;
				const int data_per_quad = node_per_face + 1;
				
				// Get sideset id
				std::string sideset_id_line = temp_data[k + 1];
				std::string sideset_id_temp;
				for (int m = 7 /*Characters in "NAME=SS" */; m < sideset_id_line.length(); ++m)
					sideset_id_temp += sideset_id_line[m];
				
				int sideset_id = 0;
				for (int m = 0; m < sideset_id_temp.length(); ++m)
					sideset_id += (sideset_id_temp[sideset_id_temp.length() - m - 1] - 48 /* ASCII TRICKS */) * pow(10.0,m);
				
				const int data_per_face = 2;
				int j = 0;
				while (( (k + 2) + j * data_per_face ) < temp_data.size() && from_string<float> (temp_float, temp_data[(k + 2) + j * data_per_face], std::dec) == true)
				{
					// Get cell to which the face belongs
					std::string face_cell_no_line = temp_data[(k + 2) + j * data_per_face];
					int face_cell_no = 0;
					for (int m = 0; m < face_cell_no_line.length(); ++m)
						face_cell_no += (face_cell_no_line[face_cell_no_line.length() - m - 1] - 48 /* ASCII TRICKS */) * pow(10.0,m);
					
					// Get ABAQUS cell face number
					std::string face_cell_face_no_line = temp_data[(k + 2) + j * data_per_face + 1];
					std::string face_cell_face_no_temp;
					for (int m = 1 /*Characters in "S" */; m < face_cell_face_no_line.length(); ++m)
						face_cell_face_no_temp += face_cell_face_no_line[m];
					
					int face_cell_face_no = 0;
					for (int m = 0; m < face_cell_face_no_temp.length(); ++m)
						face_cell_face_no += (face_cell_face_no_temp[face_cell_face_no_temp.length() - m - 1] - 48 /* ASCII TRICKS */) * pow(10.0,m);
					
					// Initilise storage variables
					std::vector <double> quad (data_per_quad);
					std::vector <double> quad_node_list (node_per_face);
					
					quad_node_list = get_global_node_numbers(face_cell_no, face_cell_face_no);
					
					// Sideset id
					quad[0] = sideset_id;
					
					// Global node numbers
					for (int m = 0; m < node_per_face; ++m)
						quad[m + 1] = quad_node_list[m];
					
					// Add to the global quad vector
					face_list.push_back(quad);
					
					++j;
				}
			}
		}
	}
	else
	{
		return false;
	}
		
	return read_in_successful;
}

//--------------------------------------------------------------------------------

std::vector <double> MeshConversion::get_global_node_numbers (const int face_cell_no, const int face_cell_face_no) 
{
	std::vector <double> quad_node_list (node_per_face);
	
	if (dimension == 2)
	{
		if (face_cell_face_no == 1)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][1];
			quad_node_list[1] = cell_list[face_cell_no - 1][2];
		}
		else if (face_cell_face_no == 2)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][2];
			quad_node_list[1] = cell_list[face_cell_no - 1][3];
		}
		else if (face_cell_face_no == 3)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][3];
			quad_node_list[1] = cell_list[face_cell_no - 1][4];
		}
		else if (face_cell_face_no == 4)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][4];
			quad_node_list[1] = cell_list[face_cell_no - 1][1];
		}
		else
		{
			std::cerr << "ERROR! Invalid face number inputted!" << std::endl;
			std::cerr << "face_cell_face_no = " << face_cell_face_no << std::endl;
		}
	}
	else // (dimension == 3)
	{
		if (face_cell_face_no == 1)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][1];
			quad_node_list[1] = cell_list[face_cell_no - 1][4];
			quad_node_list[2] = cell_list[face_cell_no - 1][3];
			quad_node_list[3] = cell_list[face_cell_no - 1][2];
		}
		else if (face_cell_face_no == 2)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][5];
			quad_node_list[1] = cell_list[face_cell_no - 1][8];
			quad_node_list[2] = cell_list[face_cell_no - 1][7];
			quad_node_list[3] = cell_list[face_cell_no - 1][6];
		}
		else if (face_cell_face_no == 3)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][1];
			quad_node_list[1] = cell_list[face_cell_no - 1][2];
			quad_node_list[2] = cell_list[face_cell_no - 1][6];
			quad_node_list[3] = cell_list[face_cell_no - 1][5];
		}
		else if (face_cell_face_no == 4)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][2];
			quad_node_list[1] = cell_list[face_cell_no - 1][3];
			quad_node_list[2] = cell_list[face_cell_no - 1][7];
			quad_node_list[3] = cell_list[face_cell_no - 1][6];
		}
		else if (face_cell_face_no == 5)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][3];
			quad_node_list[1] = cell_list[face_cell_no - 1][4];
			quad_node_list[2] = cell_list[face_cell_no - 1][8];
			quad_node_list[3] = cell_list[face_cell_no - 1][7];
		}
		else if (face_cell_face_no == 6)
		{
			quad_node_list[0] = cell_list[face_cell_no - 1][1];
			quad_node_list[1] = cell_list[face_cell_no - 1][5];
			quad_node_list[2] = cell_list[face_cell_no - 1][8];
			quad_node_list[3] = cell_list[face_cell_no - 1][4];
		}
		else
		{
			std::cerr << "ERROR! Invalid face number inputted!" << std::endl;
			std::cerr << "face_cell_face_no = " << face_cell_face_no << std::endl;
		}
	}
	
	return quad_node_list;
}

//--------------------------------------------------------------------------------

bool MeshConversion::write_out_avs_ucd (const std::string filename) 
{
	output_file_name = filename;
	bool write_out_successful = true;
	std::cout << "Writing out AVS .ucd FILE:   " << output_file_name << std::endl;
	
	std::ofstream output;
	output.open(filename.c_str());
	
	// Write out title - Note: No other commented text can be inserted below the title in a UCD file
	output << "# FEM Mesh Converter" << std::endl;
	output << "# Mesh type: AVS UCD" << std::endl;
	output << "# Input file name: " << input_file_name << std::endl;
	
	// ========================================================
	// http://people.scs.fsu.edu/~burkardt/html/ucd_format.html
	// ========================================================
	// ASCII UCD File Format
	// The input file cannot contain blank lines or lines with leading blanks. Comments, if present, must precede all data in the file. Comments within the data will cause read errors. The general order of the data is as follows:
	// 1. Numbers defining the overall structure, including the number of nodes, the number of cells, and the length of the vector of data associated with the nodes, cells, and the model.
	//     e.g. 1:
	//        <num_nodes> <num_cells> <num_ndata> <num_cdata> <num_mdata> 
	//     e.g. 2:
	//        n_elements = n_hex_cells + n_bc_quads + n_quad_cells + n_bc_edges
	//        outfile.write(str(n_nodes) + " " + str(n_elements) + " 0 0 0\n")
	// 2. For each node, its node id and the coordinates of that node in space. Node-ids must be integers, but any number including non sequential numbers can be used. Mid-edge nodes are treated like any other node.
	// 3. For each cell: its cell-id, material, cell type (hexahedral, pyramid, etc.), and the list of node-ids that correspond to each of the cell's vertices. The below table specifies the different cell types and the keyword used to represent them in the file.
	
	// Write out header
	output << node_list.size() << "\t" << (cell_list.size() + face_list.size()) << "\t0\t0\t0" << std::endl;
	
	// Write out node numbers
	for (int ii = 0; ii < node_list.size(); ++ii) // Loop over all nodes
	{
		for (int jj = 0; jj < dimension + 1; ++jj) // Loop over entries to be outputted
		{
			if (jj == 0) 	// Node number
			{
				output.precision();
				output << node_list[ii][jj] << "\t";
			}
			else 		// Node coordinates
			{
				output.width(16);
				output.setf(ios::scientific, ios::floatfield);
				output.precision(8);
				if (abs(node_list[ii][jj]) > tolerance) // invoke tolerance -> set points close to zero equal to zero
					output << static_cast<double> (node_list[ii][jj]) << "\t";
				else
					output << 0.0 << "\t";
			}
		}
		if (dimension == 2)
			output << 0.0 << "\t";
		
		output << std::endl;
		output.unsetf(ios::floatfield);
	}
	
	// Write out cell node numbers
	for (int ii = 0; ii < cell_list.size(); ++ii)
	{
		output << ii + 1 << "\t" << cell_list[ii][0] << "\t" << (dimension==2 ? "quad" : "hex") << "\t";
		for (int jj = 1; jj < node_per_cell + 1; ++jj)
			output << cell_list[ii][jj] << "\t";
		
		output << std::endl;
	}
	
	// Write out quad node numbers
	for (int ii = 0; ii < face_list.size(); ++ii)
	{
		output << ii + 1 << "\t" << face_list[ii][0] << "\t" << (dimension==2 ? "line" : "quad") << "\t";
		for (int jj = 1; jj < node_per_face + 1; ++jj)
			output << face_list[ii][jj] << "\t";
		
		output << std::endl;
	}
	
	output.close();
	
	if (write_out_successful == true)
	{
		std::cout << "Output successful!" << std::endl;
		std::cout << "   Number of nodes: " << node_list.size() << std::endl;
		std::cout << "   Number of cells: " << cell_list.size() << std::endl;
		std::cout << "   Number of boundary faces: " << face_list.size() << std::endl << std::endl;
	}
	
	return write_out_successful;
}

//--------------------------------------------------------------------------------
