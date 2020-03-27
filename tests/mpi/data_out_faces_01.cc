// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// tests DataOutFaces in parallel

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

namespace pdd
{
  // @sect3{The <code>PDDProblem</code> class template}

  template <int dim>
  class PDDProblem
  {
  public:
    PDDProblem();
    ~PDDProblem();
    void
    run();

  private:
    void
    setup_system();
    void
    output_results();

    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    // Triangulation<dim>   triangulation;
    DoFHandler<dim> dof_handler;
    FE_Q<dim>       fe;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;

    PETScWrappers::MPI::Vector locally_relevant_solution;
  };

  // The constructor and destructor of the class
  template <int dim>
  PDDProblem<dim>::PDDProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator)
    , dof_handler(triangulation)
    , fe(1)
  {}


  template <int dim>
  PDDProblem<dim>::~PDDProblem()
  {
    dof_handler.clear();
  }

  // @sect4{PDDProblem::setup_system}

  template <int dim>
  void
  PDDProblem<dim>::setup_system()
  {
    // Initialize the mesh
    std::vector<unsigned int> repetitions;
    repetitions.push_back(10);
    repetitions.push_back(12);

    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              repetitions,
                                              Point<dim>(-1., -1.),
                                              Point<dim>(1., 1.),
                                              true);
    // triangulation.refine_global (1);

    // Print out the mesh
    dof_handler.distribute_dofs(fe);

    // Initialize the solution
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    locally_relevant_solution.reinit(locally_owned_dofs, mpi_communicator);

    locally_relevant_solution = 0;

    // Apply a constant function at the boundary
    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(
      dof_handler, 0, Functions::ConstantFunction<dim>(1.0), constraints);
    constraints.close();
    constraints.distribute(locally_relevant_solution);
  }

  // Generate the outputs
  template <int dim>
  void
  PDDProblem<dim>::output_results()
  {
    // First generate an output for the cells
    const unsigned int cycle = 1;

    PETScWrappers::MPI::Vector solution(locally_owned_dofs,
                                        locally_relevant_dofs,
                                        mpi_communicator);
    solution = locally_relevant_solution;

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u");

    // make sure the following works in parallel
    data_out.build_patches();

    // on processor 0 also output some data so we can compare
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      data_out.write_vtk(deallog.get_file_stream());

    // Then generate the output for the faces.
    DataOutFaces<dim> data_out_faces;
    data_out_faces.attach_dof_handler(dof_handler);
    data_out_faces.add_data_vector(solution, "u");
    // make sure the following works in parallel
    data_out_faces.build_patches(fe.degree);

    // on processor 0 also output some data so we can compare
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      data_out_faces.write_vtk(deallog.get_file_stream());
  }
  // @sect4{PDDProblem::run}

  template <int dim>
  void
  PDDProblem<dim>::run()
  {
    setup_system();
    output_results();
  }
} // namespace pdd


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  try
    {
      pdd::PDDProblem<2> laplace_problem_2d;
      laplace_problem_2d.run();

      return 0;
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
