//----------------------------  data_out_02.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  data_out_02.cc  ---------------------------

// same as the the test with number _01, but check for block vectors


#include "../tests.h"
#include "data_out_common.cc"
#include <lac/sparsity_pattern.h>
#include <numerics/data_out.h>


std::string output_file_name = "data_out_02.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler,
            const Vector<double>  &v_node_x,
            const Vector<double>  &v_cell_x)
{
  BlockVector<double> v_node, v_cell;
  make_block_vector (v_node_x, v_node);
  make_block_vector (v_cell_x, v_cell);
  
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (v_node, "node_data");
  data_out.add_data_vector (v_cell, "cell_data");
  data_out.build_patches ();
  
  data_out.write_dx (deallog.get_file_stream());
  data_out.write_ucd (deallog.get_file_stream());  
  data_out.write_gmv (deallog.get_file_stream());
  data_out.write_tecplot (deallog.get_file_stream());
  data_out.write_vtk (deallog.get_file_stream());
  data_out.write_gnuplot (deallog.get_file_stream());

                                   // the following is only
                                   // implemented for 2d
  if (dim == 2)
    {
      data_out.write_povray (deallog.get_file_stream());
      data_out.write_eps (deallog.get_file_stream());
    }  
}


