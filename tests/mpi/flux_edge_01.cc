// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test problem in MGTools::make_flux_sparsity_pattern_edge

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

namespace Step39
{
  template <int dim>
  class InteriorPenaltyProblem
  {
  public:
    using CellInfo = MeshWorker::IntegrationInfo<dim>;

    InteriorPenaltyProblem(const FiniteElement<dim> &fe);

    void
    run(unsigned int n_steps);

  private:
    void
    setup_system();

    parallel::distributed::Triangulation<dim> triangulation;
    const MappingQ<dim>                       mapping;
    const FiniteElement<dim>                 &fe;
    DoFHandler<dim>                           dof_handler;

    IndexSet locally_relevant_set;

    TrilinosWrappers::SparseMatrix matrix;

    MGLevelObject<TrilinosWrappers::SparseMatrix> mg_matrix;

    MGLevelObject<TrilinosWrappers::SparseMatrix> mg_matrix_dg_down;
    MGLevelObject<TrilinosWrappers::SparseMatrix> mg_matrix_dg_up;
  };


  template <int dim>
  InteriorPenaltyProblem<dim>::InteriorPenaltyProblem(
    const FiniteElement<dim> &fe)
    : triangulation(MPI_COMM_WORLD,
                    Triangulation<dim>::limit_level_difference_at_vertices,
                    parallel::distributed::Triangulation<
                      dim>::construct_multigrid_hierarchy)
    , mapping(1)
    , fe(fe)
    , dof_handler(triangulation)
  {
    GridGenerator::hyper_cube_slit(triangulation, -1, 1);
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    locally_relevant_set = DoFTools::extract_locally_relevant_dofs(dof_handler);

    DynamicSparsityPattern c_sparsity(dof_handler.n_dofs(),
                                      dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, c_sparsity);
    matrix.reinit(dof_handler.locally_owned_dofs(),
                  c_sparsity,
                  MPI_COMM_WORLD,
                  true);

    const unsigned int n_levels = triangulation.n_global_levels();
    mg_matrix.resize(0, n_levels - 1);
    mg_matrix.clear_elements();
    mg_matrix_dg_up.resize(0, n_levels - 1);
    mg_matrix_dg_up.clear_elements();
    mg_matrix_dg_down.resize(0, n_levels - 1);
    mg_matrix_dg_down.clear_elements();

    for (unsigned int level = mg_matrix.min_level();
         level <= mg_matrix.max_level();
         ++level)
      {
        DynamicSparsityPattern c_sparsity(dof_handler.n_dofs(level));
        MGTools::make_flux_sparsity_pattern(dof_handler, c_sparsity, level);
        mg_matrix[level].reinit(dof_handler.locally_owned_mg_dofs(level),
                                dof_handler.locally_owned_mg_dofs(level),
                                c_sparsity,
                                MPI_COMM_WORLD,
                                true);

        if (level > 0)
          {
            DynamicSparsityPattern ci_sparsity;
            ci_sparsity.reinit(dof_handler.n_dofs(level - 1),
                               dof_handler.n_dofs(level));
            MGTools::make_flux_sparsity_pattern_edge(dof_handler,
                                                     ci_sparsity,
                                                     level);

            mg_matrix_dg_up[level].reinit(
              dof_handler.locally_owned_mg_dofs(level - 1),
              dof_handler.locally_owned_mg_dofs(level),
              ci_sparsity,
              MPI_COMM_WORLD,
              true);
            mg_matrix_dg_down[level].reinit(
              dof_handler.locally_owned_mg_dofs(level - 1),
              dof_handler.locally_owned_mg_dofs(level),
              ci_sparsity,
              MPI_COMM_WORLD,
              true);
          }
      }
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::run(unsigned int n_steps)
  {
    deallog << "Element: " << fe.get_name() << std::endl;
    for (unsigned int s = 0; s < n_steps; ++s)
      {
        deallog << "Step " << s << std::endl;
        if (s == 0)
          triangulation.refine_global(1);
        else
          {
            if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
              triangulation.begin_active(triangulation.n_levels() - 1)
                ->set_refine_flag();

            triangulation.execute_coarsening_and_refinement();
          }

        deallog << "Triangulation " << triangulation.n_active_cells()
                << " cells, " << triangulation.n_levels() << " levels"
                << std::endl;

        setup_system();
        deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, level dofs";
        for (unsigned int l = 0; l < triangulation.n_levels(); ++l)
          deallog << ' ' << dof_handler.n_dofs(l);
        deallog << std::endl;
      }
  }
} // namespace Step39



int
main(int argc, char *argv[])
{
  using namespace Step39;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  mpi_initlog(true);

  try
    {
      FE_DGQ<2>                 fe1(0);
      InteriorPenaltyProblem<2> test1(fe1);
      test1.run(20);
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
