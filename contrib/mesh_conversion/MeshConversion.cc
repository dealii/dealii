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

//--------------------------------------------------------------------------------

MeshConversion::MeshConversion (const unsigned int dimension, const std::string input_file, const std::string output_file)
    :
    tolerance (5e-16), // Used to offset Cubit tolerance error when outputting value close to zero
    dimension (dimension),
    input_file_name (input_file),
    output_file_name (output_file)
{
  greeting ();

  if (dimension == 2)
    {
      node_per_face = 2; 		// Line edge
      node_per_cell = 4;		// Quad face
      face_per_cell = 2 * dimension;	// Lines per quad
    }
  else
    if (dimension == 3)
      {
        node_per_face = 4; 		// Quad face
        node_per_cell = 8; 		// Hex cell
        face_per_cell = 2 * dimension; 	// Quads per hex
      }
    else
      {
        std::cerr << "ERROR: Chosen spatial dimension is invalid!" << std::endl;
      }
}

//--------------------------------------------------------------------------------

MeshConversion::~MeshConversion ()
{
}

//--------------------------------------------------------------------------------

bool MeshConversion::convert_mesh (void)
{
  if (!read_in_abaqus())
    {
      std::cerr << "ERROR: could not parse file." << std::endl;
      return false;
    }

  return write_out_avs_ucd();
}

//--------------------------------------------------------------------------------


// ========================================================
// http://www.codeguru.com/forum/showthread.php?t=231054
// ========================================================
template <class T> bool from_string (T& t, const std::string& s, std::ios_base& (*f) (std::ios_base&))
{
  std::istringstream iss (s);
  return ! (iss >> f >> t).fail();
}

//--------------------------------------------------------------------------------

void MeshConversion::greeting ()
{
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "************************************************************************************" << std::endl;
  std::cout << "***                         FEM Mesh conversion tool                             ***" << std::endl;
  std::cout << "************************************************************************************" << std::endl;
  std::cout << "*** Author: Jean-Paul Pelteret, Timo Heister                                     ***" << std::endl;
  std::cout << "*** Date:   June 2010, September 2011                                            ***" << std::endl;
  std::cout << "*** References:                                                                  ***" << std::endl;
  std::cout << "***    http://www.dealii.org/developer/doxygen/deal.II/structGeometryInfo.html   ***" << std::endl;
  std::cout << "***    http://people.scs.fsu.edu/~burkardt/data/ucd/ucd.html                     ***" << std::endl;
  std::cout << "***    http://www.cprogramming.com/tutorial/string.html                          ***" << std::endl;
  std::cout << "***    http://www.codeguru.com/forum/showthread.php?t=231054                     ***" << std::endl;
  std::cout << "************************************************************************************" << std::endl;
  std::cout << "*** FEATURES:                                                                    ***" << std::endl;
  std::cout << "*** Read-in file types:                                                          ***" << std::endl;
  std::cout << "***      - Abaqus inp (from Cubit)                                               ***" << std::endl;
  std::cout << "*** Write-out file types:                                                        ***" << std::endl;
  std::cout << "***      - AVS UCD                                                               ***" << std::endl;
  std::cout << "************************************************************************************" << std::endl;
  std::cout << std::endl;
}



bool MeshConversion::read_in_abaqus ()
{
  std::cout << "Reading in ABAQUS .inp FILE: " << input_file_name << std::endl;

  std::ifstream input_stream;

  input_stream.open (input_file_name.c_str(), std::ifstream::in);

  if (!input_stream.good())
    {
      std::cerr << "ERROR: file '" << input_file_name << "' not found." << std::endl;
      return false;
    }

  std::string line;

  std::getline (input_stream, line);

  while (!input_stream.eof())
    {
      if (0 == line.compare ("*HEADING")
          || 0 == line.compare (0, 2, "**")
          || 0 == line.compare (0, 5, "*PART"))
        {
          //skip header and comments
          while (!input_stream.eof())
            {
              std::getline (input_stream, line);
              if (line[0] == '*')
                goto cont;
            }
        }
      else
        if (0 == line.compare (0, 5, "*NODE"))
          {
            //header might be:
            //*NODE, NSET=ALLNODES
            //*NODE

            // contains lines in the form:
            //index, x, y, z
            int c = 0;
            while (!input_stream.eof())
              {
                std::getline (input_stream, line);
                if (line[0] == '*')
                  goto cont;

                std::vector <double> node (dimension + 1);

                istringstream iss (line);
                char comma;
                for (unsigned int i = 0; i < dimension + 1; ++i)
                  iss >> node[i] >> comma;

                node_list.push_back (node);

                ++c;
              }
          }
        else
          if (0 == line.compare (0, 8, "*ELEMENT"))
            {
              //different header formats:
              //*ELEMENT, TYPE=S4R, ELSET=EB<material id>
              //*ELEMENT, TYPE=C3D8R, ELSET=EB<material id>
              //*ELEMENT, TYPE=C3D8

              //elements itself (n=4 or n=8):
              //index, i[0], ..., i[n]

              int material = 0;

              //scan for material
              {
                string before_material = "ELSET=EB";
                size_t idx = line.find (before_material);
                if (idx != std::string::npos)
                  {
                    from_string (material, line.substr (idx + before_material.size()), std::dec);
                  }
              }

              std::getline (input_stream, line);
              while (!input_stream.eof())
                {
                  if (line[0] == '*')
                    goto cont;

                  istringstream iss (line);
                  char comma;

                  const int n_points = 2 << dimension;
                  std::vector <double> cell (1 + n_points);
                  for (int i = 0; i < 1 + n_points; ++i)
                    iss >> cell[i] >> comma;

                  //overwrite cell index from file by material
                  cell[0] = (double) material;

                  cell_list.push_back (cell);

                  std::getline (input_stream, line);
                }
            }
          else
            if (0 == line.compare (0, 8, "*SURFACE"))
              {
                //old format:
                //*SURFACE, NAME=SS<boundary indicator>
                //    <element index>,     S<face number>

                int b_indicator = 0;
                {
                  string before = "NAME=SS";
                  size_t idx = line.find (before);
                  if (idx != std::string::npos)
                    {
                      from_string (b_indicator, line.substr (idx + before.size()), std::dec);
                    }
                }

                std::getline (input_stream, line);
                while (!input_stream.eof())
                  {
                    if (line[0] == '*')
                      goto cont;

                    istringstream iss (line);
                    char comma;
                    char temp;
                    int el_idx;
                    int face_number;
                    string face_number_str;
                    iss >> el_idx >> comma >> temp >> face_number;

                    std::vector <double> quad_node_list;
                    quad_node_list = get_global_node_numbers (el_idx, face_number);
                    quad_node_list.insert (quad_node_list.begin(), b_indicator);

                    face_list.push_back (quad_node_list);

                    std::getline (input_stream, line);
                  }
              }
            else
              if (0 == line.compare (0, 6, "*ELSET"))
                {
                  //new style line/face indicators
                  //*ELSET, ELSET=SS1_S3
                  //<el1 idx>,<el2 idx>, ...
                  //*SURFACE, NAME=SS<boundary indicator>
                  //SS1_S3, S<face number>

                  std::vector<int> elements;
                  std::getline (input_stream, line);
                  while (!input_stream.eof())
                    {
                      if (line[0] == '*')
                        break;

                      istringstream iss (line);
                      char comma;
                      int elid;
                      while (!iss.eof())
                        {
                          iss >> elid >> comma;
                          elements.push_back (elid);
                        }

                      std::getline (input_stream, line);
                    }

                  if (0 != line.compare (0, 8, "*SURFACE"))
                    {
                      std::cout << "WARNING: expected SURFACE after ELSET." << std::endl;
                      goto cont;
                    }

                  int b_indicator = 0;
                  {
                    string before = "NAME=SS";
                    size_t idx = line.find (before);
                    if (idx != std::string::npos)
                      {
                        from_string (b_indicator, line.substr (idx + before.size()), std::dec);
                      }
                  }

                  std::getline (input_stream, line);

                  istringstream iss (line);
                  char comma;
                  std::string temp;
                  int face_number;
                  iss >> temp >> comma >> face_number;

                  for (unsigned int i = 0;i < elements.size();++i)
                    {
                      std::vector <double> quad_node_list;

                      quad_node_list = get_global_node_numbers (elements[i], face_number);

                      quad_node_list.insert (quad_node_list.begin(), b_indicator);

                      face_list.push_back (quad_node_list);
                    }
                }
              else
                if (0 == line.compare (0, 5, "*NSET"))
                  {
                    //skip nodesets, we have no use for them
                    while (!input_stream.eof())
                      {
                        std::getline (input_stream, line);
                        if (line[0] == '*')
                          goto cont;
                      }
                  }
                else
                  {
                    std::cout << "Warning, skipping line:" << line << std::endl;
                  }

      std::getline (input_stream, line);
    cont:
      (void) 0;
    }
  return true;
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
      else
        if (face_cell_face_no == 2)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][2];
            quad_node_list[1] = cell_list[face_cell_no - 1][3];
          }
        else
          if (face_cell_face_no == 3)
            {
              quad_node_list[0] = cell_list[face_cell_no - 1][3];
              quad_node_list[1] = cell_list[face_cell_no - 1][4];
            }
          else
            if (face_cell_face_no == 4)
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
      else
        if (face_cell_face_no == 2)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][5];
            quad_node_list[1] = cell_list[face_cell_no - 1][8];
            quad_node_list[2] = cell_list[face_cell_no - 1][7];
            quad_node_list[3] = cell_list[face_cell_no - 1][6];
          }
        else
          if (face_cell_face_no == 3)
            {
              quad_node_list[0] = cell_list[face_cell_no - 1][1];
              quad_node_list[1] = cell_list[face_cell_no - 1][2];
              quad_node_list[2] = cell_list[face_cell_no - 1][6];
              quad_node_list[3] = cell_list[face_cell_no - 1][5];
            }
          else
            if (face_cell_face_no == 4)
              {
                quad_node_list[0] = cell_list[face_cell_no - 1][2];
                quad_node_list[1] = cell_list[face_cell_no - 1][3];
                quad_node_list[2] = cell_list[face_cell_no - 1][7];
                quad_node_list[3] = cell_list[face_cell_no - 1][6];
              }
            else
              if (face_cell_face_no == 5)
                {
                  quad_node_list[0] = cell_list[face_cell_no - 1][3];
                  quad_node_list[1] = cell_list[face_cell_no - 1][4];
                  quad_node_list[2] = cell_list[face_cell_no - 1][8];
                  quad_node_list[3] = cell_list[face_cell_no - 1][7];
                }
              else
                if (face_cell_face_no == 6)
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

bool MeshConversion::write_out_avs_ucd (void)
{
  std::cout << "Writing out AVS .ucd FILE:   " << output_file_name << std::endl;

  std::ofstream output;
  output.open (output_file_name.c_str());

  // Write out title - Note: No other commented text can be inserted below the title in a UCD file
  output << "# FEM Mesh Converter" << std::endl;
  output << "# Mesh type: AVS UCD" << std::endl;
  output << "# Input file name: " << input_file_name << std::endl;

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
  for (unsigned int ii = 0; ii < node_list.size(); ++ii) // Loop over all nodes
    {
      for (unsigned int jj = 0; jj < dimension + 1; ++jj) // Loop over entries to be outputted
        {
          if (jj == 0) 	// Node number
            {
              output.precision();
              output << node_list[ii][jj] << "\t";
            }
          else 		// Node coordinates
            {
              output.width (16);
              output.setf (ios::scientific, ios::floatfield);
              output.precision (8);
              if (abs (node_list[ii][jj]) > tolerance) // invoke tolerance -> set points close to zero equal to zero
                output << static_cast<double> (node_list[ii][jj]) << "\t";
              else
                output << 0.0 << "\t";
            }
        }
      if (dimension == 2)
        output << 0.0 << "\t";

      output << std::endl;
      output.unsetf (ios::floatfield);
    }

  // Write out cell node numbers
  for (unsigned int ii = 0; ii < cell_list.size(); ++ii)
    {
      output << ii + 1 << "\t" << cell_list[ii][0] << "\t" << (dimension == 2 ? "quad" : "hex") << "\t";
      for (unsigned int jj = 1; jj < node_per_cell + 1; ++jj)
        output << cell_list[ii][jj] << "\t";

      output << std::endl;
    }

  // Write out quad node numbers
  for (unsigned int ii = 0; ii < face_list.size(); ++ii)
    {
      output << ii + 1 << "\t" << face_list[ii][0] << "\t" << (dimension == 2 ? "line" : "quad") << "\t";
      for (unsigned int jj = 1; jj < node_per_face + 1; ++jj)
        output << face_list[ii][jj] << "\t";

      output << std::endl;
    }

  output.close();

  std::cout << "Output successful!" << std::endl;
  std::cout << "   Number of nodes: " << node_list.size() << std::endl;
  std::cout << "   Number of cells: " << cell_list.size() << std::endl;
  std::cout << "   Number of boundary faces: " << face_list.size() << std::endl << std::endl;

  return true;
}

//--------------------------------------------------------------------------------
