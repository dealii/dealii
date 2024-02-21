// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check MGTransferMatrixFree and MGTransferPrebuilt for periodic boundary
// conditions

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

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

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tr, -1, 1, true);
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodicity_vector;
  GridTools::collect_periodic_faces(tr, 0, 1, 0, periodicity_vector);
  GridTools::collect_periodic_faces(tr, 2, 3, 1, periodicity_vector);
  tr.add_periodicity(periodicity_vector);
  tr.refine_global(2);
  deallog << "no. cells: " << tr.n_global_active_cells() << std::endl;

  DoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  mgdof.distribute_mg_dofs();

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(mgdof);

  // build reference
  MGTransferPrebuilt<LinearAlgebra::distributed::Vector<double>> transfer_ref(
    mg_constrained_dofs);
  transfer_ref.build(mgdof);
  deallog << "Transfer matrices: " << std::endl;
  transfer_ref.print_matrices(deallog.get_file_stream());
  deallog << std::endl;

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
      if (v3.l2_norm() > 1e-12)
        {
          // On level 0, we expect the matrix-based constraints to be wrong
          // because it cannot capture the periodicity connections with a
          // single cell
          v2_cpy.print(deallog.get_file_stream());
          v3.print(deallog.get_file_stream());
        }
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


int
main(int argc, char **argv)
{
  // no threading in this test...
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  mpi_initlog();

  check<2, double>(1);
}
