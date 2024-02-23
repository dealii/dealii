// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

#include "data_out_common.h"



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler,
           const Vector<double>  &v_node,
           const Vector<double>  &v_cell)
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(v_node, "node_data", DataOut<dim>::type_dof_data);
  data_out.add_data_vector(v_cell, "cell_data", DataOut<dim>::type_cell_data);
  data_out.build_patches();

  data_out.write_dx(deallog.get_file_stream());
  data_out.set_flags(DataOutBase::UcdFlags(true));
  data_out.write_ucd(deallog.get_file_stream());
  data_out.write_gmv(deallog.get_file_stream());
  data_out.write_tecplot(deallog.get_file_stream());
  data_out.write_vtk(deallog.get_file_stream());
  data_out.write_gnuplot(deallog.get_file_stream());
  data_out.write_deal_II_intermediate(deallog.get_file_stream());

  // the following is only
  // implemented for 2d
  if (dim == 2)
    {
      data_out.write_povray(deallog.get_file_stream());
      data_out.write_eps(deallog.get_file_stream());
    }
}
