//----------------------------  data_out_rotation_01.cc  ---------------------------
//    data_out_rotation_01.cc,v 1.1 2003/02/27 14:58:45 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  data_out_rotation_01.cc  ---------------------------

#include "../tests.h"
#include "data_out_common.cc"
#include <lac/sparsity_pattern.h>
#include <numerics/data_out_rotation.h>


std::string output_file_name = "data_out_rotation_01.output";



void
check_this (const DoFHandler<3>   &,
            const Vector<double>  &,
            const Vector<double>  &)
{
                                   // nothing to check in 3d
}


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler,
            const Vector<double>  &v_node,
            const Vector<double>  &v_cell)
{
  DataOutRotation<dim> data_out_rotation;
  data_out_rotation.attach_dof_handler (dof_handler);
  data_out_rotation.add_data_vector (v_node, "node_data");
  data_out_rotation.add_data_vector (v_cell, "cell_data");
  data_out_rotation.build_patches (4);
  
  data_out_rotation.write_dx (deallog.get_file_stream());
  data_out_rotation.write_ucd (deallog.get_file_stream());  
  data_out_rotation.write_gmv (deallog.get_file_stream());
  data_out_rotation.write_tecplot (deallog.get_file_stream());
  data_out_rotation.write_vtk (deallog.get_file_stream());
  data_out_rotation.write_gnuplot (deallog.get_file_stream());

                                   // following only implemented for
                                   // 1d+rotation=2d
  if (dim == 1)
    {
      data_out_rotation.write_povray (deallog.get_file_stream());
      data_out_rotation.write_eps (deallog.get_file_stream());
    }
}


