// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test checks that writing empty vectors does not access unallocated
// memory.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria(Triangulation<dim>::MeshSmoothing::none, true);
  GridGenerator::hyper_cube(tria, 1., 2.);
  tria.refine_global(1);

  FE_Q<dim> fe1(1);

  DoFHandler<dim> dof1(tria);
  dof1.distribute_dofs(fe1);

  Vector<double> v1(dof1.n_dofs());

  DataOut<dim> data_out;
  data_out.add_data_vector(dof1, v1, "linear");

  DataOutBase::DataOutFilter data_filter(
    DataOutBase::DataOutFilterFlags(true, false));

  data_out.write_filtered_data(data_filter);

  deallog << "Number of filtered nodes: " << data_filter.n_nodes() << std::endl;

  deallog << "ok" << std::endl;
}


int
main(int argc, char *argv[])
{
  initlog();
  test<2>();
  test<3>();

  return 0;
}
