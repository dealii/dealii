//----------------------------  data_out_faces_02.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  data_out_faces_02.cc  ---------------------------

// same as the the test with number _01, but check for block vectors

#include "data_out_common.cc"
#include <lac/sparsity_pattern.h>
#include <numerics/data_out_faces.h>


std::string output_file_name = "data_out_faces_02.output";



void
check_this (const DoFHandler<1>   &,
            const Vector<double>  &,
            const Vector<double>  &)
{
                                   // nothing to check in 1d
}


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler,
            const Vector<double>  &v_node_x,
            const Vector<double>  &v_cell_x)
{
  BlockVector<double> v_node, v_cell;
  make_block_vector (v_node_x, v_node);
  make_block_vector (v_cell_x, v_cell);
  
  DataOutFaces<dim> data_out_faces;
  data_out_faces.attach_dof_handler (dof_handler);
  data_out_faces.add_data_vector (v_node, "node_data");
  data_out_faces.add_data_vector (v_cell, "cell_data");
  data_out_faces.build_patches ();
  
  data_out_faces.write_dx (deallog.get_file_stream());
  data_out_faces.write_ucd (deallog.get_file_stream());  
  data_out_faces.write_gmv (deallog.get_file_stream());
  data_out_faces.write_tecplot (deallog.get_file_stream());
  data_out_faces.write_vtk (deallog.get_file_stream());
  data_out_faces.write_gnuplot (deallog.get_file_stream());

				   // povray and eps cannot presently
				   // write out face data
}


