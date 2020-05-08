// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// tests setup of Trilinos sparsity patterns when some processors do not have
// any cells.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include "../tests.h"

namespace Step22
{
  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem(const unsigned int degree);
    void
    run();

  private:
    void
    setup_dofs();

    const unsigned int degree;

    MPI_Comm mpi_communicator;

    SphericalManifold<dim>                    boundary;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;

    AffineConstraints<double> constraints;
    std::vector<IndexSet>     owned_partitioning;
    std::vector<IndexSet>     relevant_partitioning;
  };



  template <int dim>
  StokesProblem<dim>::StokesProblem (const unsigned int degree)
    :
    degree (degree),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator/*,
                   Triangulation<dim>::maximum_smoothing*/),
    fe (FE_Q<dim>(degree+1), dim,
        FE_Q<dim>(degree), 1),
    dof_handler (triangulation)
  {}



  template <int dim>
  void
  StokesProblem<dim>::setup_dofs()
  {
    dof_handler.distribute_dofs(fe);

    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];

    {
      owned_partitioning.clear();
      IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
      owned_partitioning.push_back(locally_owned_dofs.get_view(0, n_u));
      owned_partitioning.push_back(locally_owned_dofs.get_view(n_u, n_u + n_p));

      relevant_partitioning.clear();
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs(dof_handler,
                                              locally_relevant_dofs);
      relevant_partitioning.push_back(locally_relevant_dofs.get_view(0, n_u));
      relevant_partitioning.push_back(
        locally_relevant_dofs.get_view(n_u, n_u + n_p));
    }

    AffineConstraints<double> new_constraints;
    new_constraints.close();
    {
      TrilinosWrappers::BlockSparsityPattern bsp(owned_partitioning,
                                                 owned_partitioning,
                                                 relevant_partitioning,
                                                 mpi_communicator);
      DoFTools::make_sparsity_pattern(dof_handler,
                                      bsp,
                                      new_constraints,
                                      false,
                                      Utilities::MPI::this_mpi_process(
                                        mpi_communicator));

      bsp.compress();
    }
  }



  template <int dim>
  void
  StokesProblem<dim>::run()
  {
    Point<dim>   center;
    const double inner_radius = .5;
    const double outer_radius = 1.;

    GridGenerator::quarter_hyper_shell(
      triangulation, center, inner_radius, outer_radius, 0, true);

    triangulation.set_manifold(0, boundary);
    triangulation.set_manifold(1, boundary);
    setup_dofs();


    deallog << "OK" << std::endl;
  }
} // namespace Step22



int
main(int argc, char *argv[])
{
  try
    {
      using namespace Step22;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          std::ofstream logfile("output");
          deallog.attach(logfile, false);
          {
            StokesProblem<3> flow_problem(1);
            flow_problem.run();
          }
        }
      else
        {
          StokesProblem<3> flow_problem(1);
          flow_problem.run();
        }
    }
  catch (std::exception &exc)
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
