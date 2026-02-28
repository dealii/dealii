// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// similar to number _01, but check for auto merge

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

#include "data_out_common.h"



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler,
           const Vector<double> & v_node,
           const Vector<double> & v_cell)
{
  // Perform test only for scalar DoFHandler
  if (dof_handler.get_fe(0).n_components() > 1)
    {
      deallog.get_file_stream() << "skipped" << std::endl;
      return;
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  // Vector components
  for (unsigned int d = 0; d < dim; ++d)
    data_out.add_data_vector(v_node,
                             "vector_node_data",
                             DataOut<dim>::type_dof_data);
  for (unsigned int d = 0; d < dim; ++d)
    data_out.add_data_vector(v_cell,
                             "vector_cell_data",
                             DataOut<dim>::type_cell_data);

  // Tensor components
  for (unsigned int d = 0; d < dim * dim; ++d)
    data_out.add_data_vector(v_node,
                             "tensor_node_data",
                             DataOut<dim>::type_dof_data);
  for (unsigned int d = 0; d < dim * dim; ++d)
    data_out.add_data_vector(v_cell,
                             "tensor_cell_data",
                             DataOut<dim>::type_cell_data);

  data_out.build_patches();

  data_out.write_dx(deallog.get_file_stream());
  data_out.set_flags(DataOutBase::UcdFlags(true));
  data_out.write_ucd(deallog.get_file_stream());
  data_out.write_gmv(deallog.get_file_stream());
  data_out.write_tecplot(deallog.get_file_stream());
  data_out.write_gnuplot(deallog.get_file_stream());
  data_out.write_deal_II_intermediate(deallog.get_file_stream());

  // the following is only
  // implemented for 2d
  if (dim == 2)
    {
      data_out.write_povray(deallog.get_file_stream());
      data_out.write_eps(deallog.get_file_stream());
    }

  // VTK writer does not support tensorial output
  // For this reason a separate DataOut object is used
  DataOut<dim> data_out_vec_only;
  data_out_vec_only.attach_dof_handler(dof_handler);

  // Vector components
  for (unsigned int d = 0; d < dim; ++d)
    data_out_vec_only.add_data_vector(v_node,
                                      "vector_node_data",
                                      DataOut<dim>::type_dof_data);
  for (unsigned int d = 0; d < dim; ++d)
    data_out_vec_only.add_data_vector(v_cell,
                                      "vector_cell_data",
                                      DataOut<dim>::type_cell_data);

  data_out_vec_only.build_patches();
  data_out_vec_only.write_vtk(deallog.get_file_stream());
}
