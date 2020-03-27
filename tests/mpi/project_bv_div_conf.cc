// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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


// Check that project_boundary_values_div_conforming works for a
// distributed triangulation

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

namespace ResFlow
{
  template <int dim>
  class FluxBoundaryValues : public Function<dim>
  {
  public:
    FluxBoundaryValues()
      : Function<dim>(dim)
    {}

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &value) const override;

    virtual ~FluxBoundaryValues(){};
  };

  template <int dim>
  void
  FluxBoundaryValues<dim>::vector_value(const Point<dim> &p,
                                        Vector<double> &  bdry_flux) const
  {
    Assert(bdry_flux.size() == dim,
           ExcDimensionMismatch(bdry_flux.size(), dim));

    const double alpha = 0.3;
    const double beta  = 1;

    bdry_flux(0) = alpha * p[1] * p[1] / 2 + beta - alpha * p[0] * p[0] / 2;
    bdry_flux(1) = alpha * p[0] * p[1];
  }

  template <int dim>
  class ResFlowProblem
  {
  public:
    ResFlowProblem(const unsigned int degree);
    ~ResFlowProblem();

    void
    run();

  private:
    void
    make_grid();
    void
    setup_system();

    const unsigned int degree;

    MPI_Comm mpi_communicator;

    FESystem<dim>                             fe;
    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim>                           dof_handler;

    std::vector<IndexSet> owned_partitioning;
    std::vector<IndexSet> relevant_partitioning;

    AffineConstraints<double> constraints;

    ConditionalOStream pcout;
  };

  template <int dim>
  ResFlowProblem<dim>::ResFlowProblem(const unsigned int degree)
    : degree(degree)
    , mpi_communicator(MPI_COMM_WORLD)
    , fe(FE_RaviartThomas<dim>(degree), 1, FE_DGQ<dim>(degree), 1)
    , triangulation(mpi_communicator)
    , dof_handler(triangulation)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {}

  template <int dim>
  ResFlowProblem<dim>::~ResFlowProblem()
  {
    dof_handler.clear();
  }

  template <int dim>
  void
  ResFlowProblem<dim>::make_grid()
  {
    unsigned int nx = 1;
    unsigned int ny = nx, nz = nx;
    bool         colorize = true;

    std::vector<unsigned int> elem_dim;
    if (dim == 2)
      elem_dim = {nx, nz};
    else
      elem_dim = {nx, ny, nz};

    Point<dim> p1, p2;

    if (dim == 2)
      {
        p1 = Point<dim>(-1., -1.);
        p2 = Point<dim>(1., 1.);
      }
    else
      {
        p1 = Point<dim>(-1., -1., -1.);
        p2 = Point<dim>(1., 1., 1.);
      }

    GridGenerator::subdivided_hyper_rectangle(
      triangulation, elem_dim, p1, p2, colorize);
    triangulation.refine_global(2);
  } // end make_grid

  template <int dim>
  void
  ResFlowProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    std::vector<unsigned int> resflow_sub_blocks(2, 0);
    resflow_sub_blocks[1] = 1;
    // Renumber to yield block structure
    DoFRenumbering::block_wise(dof_handler);

    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, resflow_sub_blocks);

    const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];
    pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << " ("
          << n_u << '+' << n_p << ')' << std::endl;

    // We split up the IndexSet for locally owned and locally relevant DoFs
    // into two IndexSets based on how we want to create the block matrices
    // and vectors.
    owned_partitioning.resize(2);
    owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, n_u);
    owned_partitioning[1] =
      dof_handler.locally_owned_dofs().get_view(n_u, n_u + n_p);

    pcout << "Number of Owned dofs: " << dof_handler.n_locally_owned_dofs()
          << " (" << owned_partitioning[0].n_elements() << '+'
          << owned_partitioning[1].n_elements() << ')' << std::endl
          << std::endl;

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    relevant_partitioning.resize(2);
    relevant_partitioning[0] = locally_relevant_dofs.get_view(0, n_u);
    relevant_partitioning[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    pcout << "Number of Relevant dofs: " << locally_relevant_dofs.n_elements()
          << " (" << relevant_partitioning[0].n_elements() << '+'
          << relevant_partitioning[1].n_elements() << ')' << std::endl
          << std::endl;

    {
      constraints.reinit(locally_relevant_dofs);

      FEValuesExtractors::Vector velocities(0);
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);

      VectorTools::project_boundary_values_div_conforming(
        dof_handler, 0, FluxBoundaryValues<dim>(), 0, constraints);

      deallog << "Constraints" << std::endl;
      constraints.print(deallog.get_file_stream());
      constraints.close();
    }
  } // end setup_system


  template <int dim>
  void
  ResFlowProblem<dim>::run()
  {
    pcout << "Running with "
          << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    make_grid();
    pcout << "Grid made." << std::endl;
    setup_system();
    pcout << "System set up." << std::endl;
  }

} // end namespace ResFlow

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  const unsigned int dim = 2;

  try
    {
      using namespace ResFlow;

      Assert(dim == 2 || dim == 3, ExcNotImplemented());
      const unsigned int  press_order = 1;
      ResFlowProblem<dim> resflow(press_order);
      resflow.run();
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

  deallog << "OK" << std::endl;
  return 0;
}
