// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// output the vertex numbering in a vtk file

#include "../tests.h"

// all include files needed for the program

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <string>

template <int dim, int spacedim>
void
test(std::string filename)
{
  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim>        gi;

  gi.attach_triangulation(triangulation);
  std::ifstream in(filename);
  gi.read_ucd(in);

  FE_Q<dim, spacedim>       fe(1);
  DoFHandler<dim, spacedim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  // Output the vertex numbering
  Vector<double> numbering(dof_handler.n_dofs());
  for (unsigned int i = 0; i < numbering.size(); ++i)
    numbering(i) = i;

  DataOut<dim, spacedim> dataout;
  dataout.add_data_vector(dof_handler, numbering, "numbering");
  dataout.build_patches();
  dataout.write_vtk(deallog.get_file_stream());
}



int
main()
{
  initlog();
  deallog << "Test<1,2>" << std::endl;
  test<1, 2>(SOURCE_DIR "/grids/circle_2.inp");

  deallog << "Test<2,3>" << std::endl;
  test<2, 3>(SOURCE_DIR "/grids/sphere_2.inp");

  return 0;
}
