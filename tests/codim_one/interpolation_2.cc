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



// calculates the surface of a sphere

#include "../tests.h"

// all include files needed for the program

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
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

  // Now we interpolate the constant function on the mesh, and check
  // that this is consistent with what we expect.
  Vector<double> interpolated_one(dof_handler.n_dofs());

  Functions::CosineFunction<spacedim> cosine;
  VectorTools::interpolate(dof_handler, cosine, interpolated_one);

  DataOut<dim, spacedim> dataout;
  dataout.add_data_vector(dof_handler, interpolated_one, "numbering");
  dataout.build_patches();
  dataout.write_vtk(deallog.get_file_stream());
}



int
main()
{
  initlog();

  deallog << "Test<1,2>" << std::endl;
  test<1, 2>(SOURCE_DIR "/grids/circle_2.inp");

  deallog << "Test<1,2>" << std::endl;
  test<2, 3>(SOURCE_DIR "/grids/sphere_2.inp");

  return 0;
}
