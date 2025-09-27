// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// DataOut::build_patches appeared to have a problem with outputting
// lines in 2d where nodes were numbered differently when writing data
// vectors as opposed to writing node locations. in the end this
// turned out to be a feature: the mesh was a circle of lines, so
// there are equally many cells as their were nodes, and consequently
// DataOut assumed that it had cell_data, rather than
// dof_data. passing the correct argument fixed the problem, but it
// won't hurt to have this test anyway.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int dim = 1;

  Triangulation<dim, dim + 1> triangulation;
  FE_Q<dim, dim + 1>          fe(1);
  DoFHandler<dim, dim + 1>    dof_handler(triangulation);
  Vector<double>              soln;

  GridIn<dim, dim + 1> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream fname(SOURCE_DIR "/grids/square.msh");
  grid_in.read_msh(fname);

  dof_handler.distribute_dofs(fe);
  soln.reinit(dof_handler.n_dofs());
  soln = 0;
  for (unsigned int i = 0; i < soln.size(); ++i)
    soln(i) = i;
  DataOut<dim, dim + 1> data_out;
  data_out.attach_dof_handler(dof_handler);

  data_out.add_data_vector(soln,
                           "scalar_data",
                           DataOut<dim, dim + 1>::type_dof_data);
  data_out.build_patches();
  data_out.write_vtk(deallog.get_file_stream());

  return 0;
}
