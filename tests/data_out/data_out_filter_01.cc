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

// This test shows the ability of DataOutFilter to merge vertex duplicates
// taking into account floating point inprecision. The test creates a distorted
// triangulation with vertices very close to each other that are subsequently
// merged by the DataOutFilter class. Note that while the test case will likely
// never happen, there are cases in which vertices during mesh output have
// slightly different positions that would not be merged if the points would
// simply be compared for equality.

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
  for (unsigned int i = 0; i < v1.size(); ++i)
    v1(i) = i;

  DataOut<dim> data_out;
  data_out.add_data_vector(dof1, v1, "linear");
  data_out.build_patches();

  DataOutBase::DataOutFilter data_filter(
    DataOutBase::DataOutFilterFlags(true, false));

  data_out.write_filtered_data(data_filter);

  deallog << "Number of filtered nodes: " << data_filter.n_nodes() << std::endl;

  DataOutBase::DataOutFilter data_filter2(
    DataOutBase::DataOutFilterFlags(false, false));
  data_out.write_filtered_data(data_filter2);

  deallog << "Number of unfiltered nodes: " << data_filter2.n_nodes()
          << std::endl;

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
