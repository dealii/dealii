// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// same as data_out_faces_01 but without attaching a dof handler

#include "../tests.h"
#include "data_out_common.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out_faces.h>


std::string output_file_name = "output";



void
my_check_this (const DoFHandler<1> &,
               const Vector<double> &,
               const Vector<double> &)
{
  // nothing to check in 1d
}


template <int dim>
void
my_check_this (const DoFHandler<dim> &dof_handler,
               const Vector<double>  &v_node,
               const Vector<double>  &v_cell)
{
  DataOutFaces<dim> data_out_faces;
  data_out_faces.add_data_vector (dof_handler, v_node, "node_data");
  data_out_faces.add_data_vector (v_cell, "cell_data");
  data_out_faces.build_patches ();

  data_out_faces.write_dx (deallog.get_file_stream());
  data_out_faces.set_flags (DataOutBase::UcdFlags(true));
  data_out_faces.write_ucd (deallog.get_file_stream());
  data_out_faces.write_gmv (deallog.get_file_stream());
  data_out_faces.write_tecplot (deallog.get_file_stream());
  data_out_faces.write_vtk (deallog.get_file_stream());
  data_out_faces.write_gnuplot (deallog.get_file_stream());
  data_out_faces.write_deal_II_intermediate (deallog.get_file_stream());

  // povray and eps cannot presently
  // write out face data
}


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler,
            const Vector<double>  &v_node,
            const Vector<double>  &v_cell)
{
  // since we can't forward declare
  // check_this in this file (it is forward
  // declared in data_out_common.h), we
  // also can't make the driver file aware of
  // the overload for 1d. to avoid linker
  // errors, we can consequently not overload
  // check_this, and need this forwarder
  // function
  my_check_this (dof_handler, v_node, v_cell);
}
