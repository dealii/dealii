// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check MGTransferMatrixFree by comparison with MGTransferPrebuilt on a
// series of meshes with uniform meshes for FE_Q

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"


template <int dim, typename Number>
void
check(const unsigned int fe_degree)
{
  FE_Q<dim> fe(fe_degree);
  deallog << "FE: " << fe.get_name() << std::endl;

  // run a few different sizes...
  unsigned int sizes[] = {1, 2, 3, 4, 5, 6, 8};
  for (unsigned int cycle = 0; cycle < sizeof(sizes) / sizeof(unsigned int);
       ++cycle)
    {
      unsigned int n_refinements = 0;
      unsigned int n_subdiv      = sizes[cycle];
      if (n_subdiv > 1)
        while (n_subdiv % 2 == 0)
          {
            n_refinements += 1;
            n_subdiv /= 2;
          }
      n_refinements += 3 - dim;
      if (fe_degree < 3)
        n_refinements += 1;
      parallel::distributed::Triangulation<dim> tr(
        MPI_COMM_WORLD,
        Triangulation<dim>::limit_level_difference_at_vertices,
        parallel::distributed::Triangulation<
          dim>::construct_multigrid_hierarchy);
      GridGenerator::subdivided_hyper_cube(tr, n_subdiv);
      tr.refine_global(n_refinements);
      deallog << "no. cells: " << tr.n_global_active_cells() << std::endl;

      DoFHandler<dim> mgdof(tr);
      mgdof.distribute_dofs(fe);
      mgdof.distribute_mg_dofs();

      MGConstrainedDoFs mg_constrained_dofs;
      mg_constrained_dofs.initialize(mgdof);
      mg_constrained_dofs.make_zero_boundary_constraints(mgdof, {0});

      // build reference
      MGTransferPrebuilt<LinearAlgebra::distributed::Vector<double>>
        transfer_ref(mg_constrained_dofs);
      transfer_ref.build(mgdof);

      // build matrix-free transfer
      MGTransferMF<dim, Number> transfer(mg_constrained_dofs);
      transfer.build(mgdof);

      // check prolongation for all levels using random vector
      for (unsigned int level = 1;
           level < mgdof.get_triangulation().n_global_levels();
           ++level)
        {
          LinearAlgebra::distributed::Vector<Number> v1, v2;
          LinearAlgebra::distributed::Vector<double> v1_cpy, v2_cpy, v3;
          v1.reinit(mgdof.locally_owned_mg_dofs(level - 1), MPI_COMM_WORLD);
          v2.reinit(mgdof.locally_owned_mg_dofs(level), MPI_COMM_WORLD);
          v3.reinit(mgdof.locally_owned_mg_dofs(level), MPI_COMM_WORLD);
          for (unsigned int i = 0; i < v1.locally_owned_size(); ++i)
            v1.local_element(i) = random_value<double>();
          v1_cpy = v1;
          transfer.prolongate(level, v2, v1);
          transfer_ref.prolongate(level, v3, v1_cpy);
          v2_cpy = v2;
          v3 -= v2_cpy;
          deallog << "Diff prolongate   l" << level << ": " << v3.l2_norm()
                  << std::endl;
        }

      // check restriction for all levels using random vector
      for (unsigned int level = 1;
           level < mgdof.get_triangulation().n_global_levels();
           ++level)
        {
          LinearAlgebra::distributed::Vector<Number> v1, v2;
          LinearAlgebra::distributed::Vector<double> v1_cpy, v2_cpy, v3;
          v1.reinit(mgdof.locally_owned_mg_dofs(level), MPI_COMM_WORLD);
          v2.reinit(mgdof.locally_owned_mg_dofs(level - 1), MPI_COMM_WORLD);
          v3.reinit(mgdof.locally_owned_mg_dofs(level - 1), MPI_COMM_WORLD);
          for (unsigned int i = 0; i < v1.locally_owned_size(); ++i)
            v1.local_element(i) = random_value<double>();
          v1_cpy = v1;
          transfer.restrict_and_add(level, v2, v1);
          transfer_ref.restrict_and_add(level, v3, v1_cpy);
          v2_cpy = v2;
          v3 -= v2_cpy;
          deallog << "Diff restrict     l" << level << ": " << v3.l2_norm()
                  << std::endl;

          v2 = 1.;
          v3 = 1.;
          transfer.restrict_and_add(level, v2, v1);
          transfer_ref.restrict_and_add(level, v3, v1_cpy);
          v2_cpy = v2;
          v3 -= v2_cpy;
          deallog << "Diff restrict add l" << level << ": " << v3.l2_norm()
                  << std::endl;
        }
      deallog << std::endl;
    }
}


int
main(int argc, char **argv)
{
  // no threading in this test...
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  mpi_initlog();

  check<2, double>(1);
  check<2, double>(3);
  check<3, double>(1);
  check<3, double>(3);
  check<2, float>(2);
  check<3, float>(2);
}
