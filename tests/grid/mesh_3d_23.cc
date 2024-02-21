// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like mesh_3d_22, but further reduced: when creating output with
// DataOut using a MappingQ(3) on a mesh with flipped cells, we get
// bogus output at the interior face.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"


int
main()
{
  initlog();

  Triangulation<3> triangulation;
  GridIn<3>        grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream inputStream(SOURCE_DIR "/grids/mesh.msh");
  grid_in.read_msh(inputStream);

  MappingQ<3>   mapping(3);
  FE_Q<3>       fe(3);
  DoFHandler<3> dofh(triangulation);

  dofh.distribute_dofs(fe);

  Vector<double> x(dofh.n_dofs());
  DataOut<3>     data_out;

  data_out.attach_dof_handler(dofh);
  data_out.add_data_vector(x, "u");
  data_out.build_patches(mapping, 3);

  data_out.write_gnuplot(deallog.get_file_stream());
}
