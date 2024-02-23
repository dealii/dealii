// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this test used to fail somewhere in the 3d parts of MappingQ for fourth
// order mappings. it later turned out that the weird error was simply that
// something wasn't implemented, which is now the case

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



int
main()
{
  initlog();

  Triangulation<3> coarse_grid;
  GridGenerator::hyper_cube(coarse_grid, -1, 1);
  coarse_grid.refine_global(1);

  FE_Q<3>       fe(1);
  DoFHandler<3> dof_handler(coarse_grid);
  dof_handler.distribute_dofs(fe);

  MappingQ<3>     mapping(4);
  QGauss<2>       q_face(3);
  FEFaceValues<3> fe_face_values(mapping,
                                 fe,
                                 q_face,
                                 update_normal_vectors | update_JxW_values);
  fe_face_values.reinit(dof_handler.begin_active(), 0);

  // if we got here, we got past the previous
  // abort
  deallog << "OK" << std::endl;
}
