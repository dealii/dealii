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


// check MGTransferBlockMatrixFree::copy_[from|to]_mg by comparison to non-block
// MGTransferMatrixFree. This is similar to transfer_05 but using
// a separate DoFHandler for each block

#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <algorithm>

#include "../tests.h"



template <int dim, typename Number>
void
check(const unsigned int fe_degree)
{
  const unsigned int nb = 2;

  FE_Q<dim> fe_1(fe_degree + 1);
  FE_Q<dim> fe_2(fe_degree);

  deallog << "FE: " << fe_1.get_name() << std::endl;
  deallog << "FE: " << fe_2.get_name() << std::endl;

  // run a few different sizes...
  unsigned int sizes[] = {1, 2, 3};
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

      // adaptive refinement into a circle
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tr.begin_active();
           cell != tr.end();
           ++cell)
        if (cell->is_locally_owned() && cell->center().norm() < 0.5)
          cell->set_refine_flag();
      tr.execute_coarsening_and_refinement();
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tr.begin_active();
           cell != tr.end();
           ++cell)
        if (cell->is_locally_owned() && cell->center().norm() > 0.3 &&
            cell->center().norm() < 0.4)
          cell->set_refine_flag();
      tr.execute_coarsening_and_refinement();
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tr.begin_active();
           cell != tr.end();
           ++cell)
        if (cell->is_locally_owned() && cell->center().norm() > 0.33 &&
            cell->center().norm() < 0.37)
          cell->set_refine_flag();
      tr.execute_coarsening_and_refinement();

      deallog << "no. cells: " << tr.n_global_active_cells() << std::endl;

      DoFHandler<dim> mgdof_1(tr);
      mgdof_1.distribute_dofs(fe_1);
      mgdof_1.distribute_mg_dofs();

      DoFHandler<dim> mgdof_2(tr);
      mgdof_2.distribute_dofs(fe_2);
      mgdof_2.distribute_mg_dofs();

      const std::vector<const DoFHandler<dim> *> mgdof_ptr{&mgdof_1, &mgdof_2};

      std::vector<MGConstrainedDoFs> mg_constrained_dofs_vector(2);
      for (unsigned int i = 0; i < mgdof_ptr.size(); ++i)
        {
          mg_constrained_dofs_vector[i].initialize(*mgdof_ptr[i]);
          mg_constrained_dofs_vector[i].make_zero_boundary_constraints(
            *mgdof_ptr[i], {0});
        }

      // build reference non-block
      std::vector<MGTransferMatrixFree<dim, Number>> transfer_ref;
      for (unsigned int b = 0; b < nb; ++b)
        {
          transfer_ref.emplace_back(mg_constrained_dofs_vector[b]);
          transfer_ref[b].build(*mgdof_ptr[b]);
        }

      // build matrix-free block transfer
      MGTransferBlockMatrixFree<dim, Number> transfer(
        mg_constrained_dofs_vector);
      transfer.build(mgdof_ptr);

      MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> lbv(
        0, tr.n_global_levels() - 1);
      LinearAlgebra::distributed::BlockVector<Number> bv(nb);

      MGLevelObject<LinearAlgebra::distributed::Vector<Number>> lv(
        0, tr.n_global_levels() - 1);
      LinearAlgebra::distributed::Vector<Number> v(nb);

      // initialize
      for (unsigned int b = 0; b < nb; ++b)
        bv.block(b).reinit(mgdof_ptr[b]->locally_owned_dofs(), MPI_COMM_WORLD);

      for (unsigned int l = lbv.min_level(); l <= lbv.max_level(); ++l)
        {
          lbv[l].reinit(nb);
          for (unsigned int b = 0; b < nb; ++b)
            lbv[l].block(b).reinit(mgdof_ptr[b]->locally_owned_mg_dofs(l),
                                   MPI_COMM_WORLD);

          lbv[l].collect_sizes();

          // set values:
          for (unsigned int b = 0; b < nb; ++b)
            for (unsigned int i = 0; i < lbv[l].block(b).locally_owned_size();
                 ++i)
              lbv[l].block(b).local_element(i) = random_value<double>();

          lbv[l].compress(VectorOperation::insert);
        }

      // check copy_from_mg
      transfer.copy_from_mg(mgdof_ptr, bv, lbv);
      for (unsigned int b = 0; b < nb; ++b)
        {
          v.reinit(mgdof_ptr[b]->locally_owned_dofs(), MPI_COMM_WORLD);

          for (unsigned int l = lv.min_level(); l <= lv.max_level(); ++l)
            {
              lv[l].reinit(mgdof_ptr[b]->locally_owned_mg_dofs(l),
                           MPI_COMM_WORLD);
              lv[l] = lbv[l].block(b);
            }

          transfer_ref[b].copy_from_mg(*mgdof_ptr[b], v, lv);
          v -= bv.block(b);
          deallog << "Diff copy_from_mg b" << b << ": " << v.l2_norm()
                  << std::endl;
        }

      // check copy_to_mg
      // use un-initialized level-block vector to make sure that it will be
      // set correctly in copy_to_mg().
      MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> lbv2(
        0, tr.n_global_levels() - 1);
      for (unsigned int b = 0; b < nb; ++b)
        for (unsigned int i = 0; i < bv.block(b).locally_owned_size(); ++i)
          bv.block(b).local_element(i) = random_value<double>();

      transfer.copy_to_mg(mgdof_ptr, lbv2, bv);
      // Also check that the block vector has its (global) size set on each
      // level:
      for (unsigned int l = lv.min_level(); l <= lv.max_level(); ++l)
        {
          unsigned int total_size = 0;
          for (unsigned int b = 0; b < nb; ++b)
            total_size += lbv2[l].block(b).size();

          AssertThrow(total_size == lbv2[l].size(),
                      ExcDimensionMismatch(total_size, lbv2[l].size()));
        }

      // Finally check the difference:
      for (unsigned int b = 0; b < nb; ++b)
        {
          v.reinit(mgdof_ptr[b]->locally_owned_dofs(), MPI_COMM_WORLD);

          v = bv.block(b);
          transfer_ref[b].copy_to_mg(*mgdof_ptr[b], lv, v);
          for (unsigned int l = lv.min_level(); l <= lv.max_level(); ++l)
            {
              lv[l] -= lbv2[l].block(b);
              deallog << "Diff copy_to_mg   l" << l << ": " << lv[l].l2_norm()
                      << std::endl;
            }
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
