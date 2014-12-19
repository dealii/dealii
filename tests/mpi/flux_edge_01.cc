// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test problem in MGTools::make_flux_sparsity_pattern_edge

#include "../tests.h"
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <iostream>
#include <fstream>

namespace Step39
{
  using namespace dealii;

  template <int dim>
  class InteriorPenaltyProblem
  {
  public:
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;

    InteriorPenaltyProblem(const FiniteElement<dim> &fe);

    void run(unsigned int n_steps);

  private:
    void setup_system ();

    parallel::distributed::Triangulation<dim>        triangulation;
    const MappingQ1<dim>      mapping;
    const FiniteElement<dim> &fe;
    DoFHandler<dim>           dof_handler;

    IndexSet locally_relevant_set;

    TrilinosWrappers::SparseMatrix matrix;

    MGLevelObject<TrilinosWrappers::SparseMatrix> mg_matrix;

    MGLevelObject<TrilinosWrappers::SparseMatrix> mg_matrix_dg_down;
    MGLevelObject<TrilinosWrappers::SparseMatrix> mg_matrix_dg_up;
  };


  template <int dim>
  InteriorPenaltyProblem<dim>::InteriorPenaltyProblem(const FiniteElement<dim> &fe)
    :
    triangulation (MPI_COMM_WORLD,Triangulation<dim>::
                   limit_level_difference_at_vertices,
                   parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy),
    mapping(),
    fe(fe),
    dof_handler(triangulation)
  {
    GridGenerator::hyper_cube_slit(triangulation, -1, 1);
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs (fe);

    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_set);

    CompressedSimpleSparsityPattern c_sparsity(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, c_sparsity);
    matrix.reinit(dof_handler.locally_owned_dofs(), c_sparsity, MPI_COMM_WORLD, true);

    const unsigned int n_levels = triangulation.n_global_levels();
    mg_matrix.resize(0, n_levels-1);
    mg_matrix.clear();
    mg_matrix_dg_up.resize(0, n_levels-1);
    mg_matrix_dg_up.clear();
    mg_matrix_dg_down.resize(0, n_levels-1);
    mg_matrix_dg_down.clear();

    for (unsigned int level=mg_matrix.min_level();
         level<=mg_matrix.max_level(); ++level)
      {
        CompressedSimpleSparsityPattern c_sparsity(dof_handler.n_dofs(level));
        MGTools::make_flux_sparsity_pattern(dof_handler, c_sparsity, level);
        mg_matrix[level].reinit(dof_handler.locally_owned_mg_dofs(level),
                                dof_handler.locally_owned_mg_dofs(level),
                                c_sparsity,
                                MPI_COMM_WORLD, true);

        if (level>0)
          {
            CompressedSimpleSparsityPattern ci_sparsity;
            ci_sparsity.reinit(dof_handler.n_dofs(level-1), dof_handler.n_dofs(level));
            MGTools::make_flux_sparsity_pattern_edge(dof_handler, ci_sparsity, level);

            mg_matrix_dg_up[level].reinit(dof_handler.locally_owned_mg_dofs(level-1),
                                          dof_handler.locally_owned_mg_dofs(level),
                                          ci_sparsity,
                                          MPI_COMM_WORLD, true);
            mg_matrix_dg_down[level].reinit(dof_handler.locally_owned_mg_dofs(level-1),
                                            dof_handler.locally_owned_mg_dofs(level),
                                            ci_sparsity,
                                            MPI_COMM_WORLD, true);
          }
      }
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::run(unsigned int n_steps)
  {
    deallog << "Element: " << fe.get_name() << std::endl;
    for (unsigned int s=0; s<n_steps; ++s)
      {
        deallog << "Step " << s << std::endl;
        if (s==0)
          triangulation.refine_global(1);
        else
          {
            if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
              triangulation.begin_active(triangulation.n_levels()-1)->set_refine_flag();

            triangulation.execute_coarsening_and_refinement ();
          }

        deallog << "Triangulation "
                << triangulation.n_active_cells() << " cells, "
                << triangulation.n_levels() << " levels" << std::endl;

        setup_system();
        deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, level dofs";
        for (unsigned int l=0; l<triangulation.n_levels(); ++l)
          deallog << ' ' << dof_handler.n_dofs(l);
        deallog << std::endl;
      }
  }
}



int main(int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  mpi_initlog(true);

  using namespace dealii;
  using namespace Step39;

  try
    {
      FE_DGQ<2> fe1(0);
      InteriorPenaltyProblem<2> test1(fe1);
      test1.run(20);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
