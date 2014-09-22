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


#include "../tests.h"
#include "data_out_common.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out_stack.h>


std::string output_file_name = "output";



template <int dim>
void
check_this (const DoFHandler<dim> &,
            const Vector<double>  &,
            const Vector<double>  &)
{
  // 3d would generate 4d data, which
  // we don't presently support
  //
  // output for 2d+time is not
  // presently implemented
}



template <>
void
check_this<1> (const DoFHandler<1> &dof_handler,
               const Vector<double>  &v_node,
               const Vector<double>  &v_cell)
{
  const unsigned int dim = 1;
  
  DataOutStack<dim> data_out_stack;
  data_out_stack.declare_data_vector ("node_data",
                                      DataOutStack<dim>::dof_vector);
  data_out_stack.declare_data_vector ("cell_data",
                                      DataOutStack<dim>::cell_vector);


  data_out_stack.new_parameter_value (1., 1.);
  data_out_stack.attach_dof_handler (dof_handler);
  data_out_stack.add_data_vector (v_node, "node_data");
  data_out_stack.add_data_vector (v_cell, "cell_data");
  data_out_stack.build_patches ();
  data_out_stack.finish_parameter_value ();

  Vector<double> vn1 (v_node.size());
  vn1 = v_node;
  vn1 *= 2.;
  Vector<double> vc1 (v_cell.size());
  vc1 = v_cell;
  vc1 *= 3.;
  data_out_stack.new_parameter_value (1., 1.);
  data_out_stack.attach_dof_handler (dof_handler);
  data_out_stack.add_data_vector (vn1, "node_data");
  data_out_stack.add_data_vector (vc1, "cell_data");
  data_out_stack.build_patches ();
  data_out_stack.finish_parameter_value ();

  data_out_stack.write_dx (deallog.get_file_stream());
  data_out_stack.set_flags (DataOutBase::UcdFlags(true));
  data_out_stack.write_ucd (deallog.get_file_stream());
  data_out_stack.write_gmv (deallog.get_file_stream());
  data_out_stack.write_tecplot (deallog.get_file_stream());
  data_out_stack.write_vtk (deallog.get_file_stream());
  data_out_stack.write_gnuplot (deallog.get_file_stream());
  data_out_stack.write_deal_II_intermediate (deallog.get_file_stream());

  // the following is only
  // implemented for 2d (=1d+time)
  if (dim == 1)
    {
      data_out_stack.write_povray (deallog.get_file_stream());
      data_out_stack.write_eps (deallog.get_file_stream());
    }
}


