// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test for interpolate_to_different_mesh with dealii::DoFHandler
// checks that interpolation from coarse to fine and back to coarse
// is the identity operation

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include "../tests.h"



template <unsigned int spacedim>
void
check(const unsigned int refinement_1, const unsigned int refinement_2)
{
  MappingQ<spacedim> mapping(1);

  Triangulation<spacedim> tria_1, tria_2;
  GridGenerator::hyper_cube(tria_1);
  GridGenerator::hyper_cube(tria_2);

  tria_1.refine_global(refinement_1);
  tria_2.refine_global(refinement_2);

  DoFHandler<spacedim> dof_handler_1(tria_1);
  DoFHandler<spacedim> dof_handler_2(tria_2);

  FE_Q<spacedim> fe(1);
  dof_handler_1.distribute_dofs(fe);
  dof_handler_2.distribute_dofs(fe);

  Vector<double> u_1(dof_handler_1.n_dofs()), u_2(dof_handler_2.n_dofs()),
    u_1_back(dof_handler_1.n_dofs()), d(dof_handler_1.n_dofs());
  for (unsigned int i1 = 0; i1 < u_1.size(); ++i1)
    u_1[i1] = 2. * rand() / RAND_MAX - 1.;

  VectorTools::interpolate_to_different_mesh(dof_handler_1,
                                             u_1,
                                             dof_handler_2,
                                             u_2);

  VectorTools::interpolate_to_different_mesh(dof_handler_2,
                                             u_2,
                                             dof_handler_1,
                                             u_1_back);

  d = u_1;
  d.sadd(-1.0, u_1_back);
  deallog << "distance=" << d.l2_norm() << std::endl;
}

int
main(int argc, char **argv)
{
  initlog();

  srand(static_cast<unsigned int>(clock()));

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog
    << "### 2D-Case, first cell unrefined, second cell refined once###\n\n";
  check<2>(0, 1);
  deallog << std::endl;

  deallog
    << "### 2D-Case, first cell refined once, second cell refined three times###\n\n";
  check<2>(1, 3);
  deallog << std::endl;

  deallog
    << "### 3D-Case, first cell unrefined, second cell refined once###\n\n";
  check<3>(0, 1);
  deallog << std::endl;

  deallog
    << "### 3D-Case, first cell refined once, second cell refined three times###\n\n";
  check<3>(1, 3);
  deallog << std::endl;
}
