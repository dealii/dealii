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


// same as data_out_rotation_01.cc, but without attaching a dof handler

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out_rotation.h>

#include "../tests.h"

#include "data_out_common.h"



void
my_check_this(const DoFHandler<3> &,
              const Vector<double> &,
              const Vector<double> &)
{
  // nothing to check in 3d
}


template <int dim>
void
my_check_this(const DoFHandler<dim> &dof_handler,
              const Vector<double>  &v_node,
              const Vector<double>  &v_cell)
{
  DataOutRotation<dim> data_out_rotation;
  data_out_rotation.add_data_vector(dof_handler, v_node, "node_data");
  data_out_rotation.add_data_vector(v_cell, "cell_data");
  data_out_rotation.build_patches(4);

  data_out_rotation.write_dx(deallog.get_file_stream());
  data_out_rotation.set_flags(DataOutBase::UcdFlags(true));
  data_out_rotation.write_ucd(deallog.get_file_stream());
  data_out_rotation.write_gmv(deallog.get_file_stream());
  data_out_rotation.write_tecplot(deallog.get_file_stream());
  data_out_rotation.write_vtk(deallog.get_file_stream());
  data_out_rotation.write_gnuplot(deallog.get_file_stream());
  data_out_rotation.write_deal_II_intermediate(deallog.get_file_stream());

  // following only implemented for 1d+rotation=2d
  if (dim == 1)
    {
      data_out_rotation.write_povray(deallog.get_file_stream());
      data_out_rotation.write_eps(deallog.get_file_stream());
    }
}


template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler,
           const Vector<double>  &v_node,
           const Vector<double>  &v_cell)
{
  // since we can't forward declare check_this in this file (it is forward
  // declared in data_out_common.h), we also can't make the driver file aware
  // of the overload for 1d. to avoid linker errors, we can consequently not
  // overload check_this, and need this forwarder function
  my_check_this(dof_handler, v_node, v_cell);
}
