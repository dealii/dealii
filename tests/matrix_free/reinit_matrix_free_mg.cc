// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Verify that MatrixFree::reinit works without errors also when we have
// processes without any cells on a larger mesh (we need 17 processes to have
// a 8x8 mesh with one process not having a cell on level 1).


#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <fstream>
#include <iostream>

#include "../tests.h"



template <int dim>
void
do_test(const unsigned int n_refinements)
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);

  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  deallog << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

  AffineConstraints<double> constraints;
  constraints.close();

  MappingQ1<dim> mapping;

  for (unsigned int level = 0; level <= tria.n_global_levels(); ++level)
    {
      typename MatrixFree<dim, double>::AdditionalData data;
      data.initialize_mapping = false;
      if (level < tria.n_global_levels())
        data.mg_level = level;

      MatrixFree<dim, double> matrix_free;
      matrix_free.reinit(mapping, dof_handler, constraints, QGauss<1>(2), data);

      deallog << "Level " << data.mg_level << " OK" << std::endl;
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  mpi_initlog();

  do_test<2>(2);
  do_test<2>(3);
  do_test<2>(4);
}
