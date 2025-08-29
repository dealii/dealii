// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found incd
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Step-40 on simplex meshes, reading the grid from separate partitioned GMSH
// files.
//
// The key modifications with respect to the original Step-40 tutorial are:
//
// - Grid generation is performed using the Gmsh C++ API.
// - The mesh files are created using Gmsh’s `.geo` description and exported
//   as `.msh` files containing triangular (simplex) elements.
//
// As in the quadrilateral case, partitioning is performed using Gmsh's
// built-in command-line utilities or C++ API with flags:
//
//     fname.msh -part N -part_ghosts -part_split -save
//
// Each MPI process reads only its corresponding partitioned simplex mesh file.
// Ghost cells are generated automatically during partitioning, and no global
// coarse mesh is required to exist on all processes.
//
// This setup ensures scalable parallel performance and efficient memory use
// for simulations involving unstructured triangular elements.


#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA

namespace Step40
{
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem();
    void
    run();

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    output_results() const;

    MPI_Comm                                       mpi_communicator;
    MappingFE<dim>                                 mapping;
    parallel::fullydistributed::Triangulation<dim> triangulation;
    FE_SimplexP<dim>                               fe;
    DoFHandler<dim>                                dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector       locally_relevant_solution;
    LA::MPI::Vector       system_rhs;
  };

  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , mapping(FE_SimplexP<dim>(1))
    , triangulation(mpi_communicator)
    , fe(2)
    , dof_handler(triangulation)
  {}

  template <int dim>
  void
  LaplaceProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    MPI_Barrier(mpi_communicator);

    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);

    constraints.clear();
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(
      mapping, dof_handler, 11, Functions::ZeroFunction<dim>(), constraints);
    constraints.close();

    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.locally_owned_dofs(),
                                               mpi_communicator,
                                               locally_relevant_dofs);
    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);
  }

  template <int dim>
  void
  LaplaceProblem<dim>::assemble_system()
  {
    const QGaussSimplex<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(mapping,
                            fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      {
        cell_matrix = 0.;
        cell_rhs    = 0.;
        fe_values.reinit(cell);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double rhs_value =
              (fe_values.quadrature_point(q_point)[1] >
                   0.5 +
                     0.25 * std::sin(4.0 * numbers::PI *
                                     fe_values.quadrature_point(q_point)[0]) ?
                 1. :
                 -1.);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
                                       fe_values.shape_grad(j, q_point) *
                                       fe_values.JxW(q_point);

                cell_rhs(i) += rhs_value *                         //
                               fe_values.shape_value(i, q_point) * //
                               fe_values.JxW(q_point);
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }

  template <int dim>
  void
  LaplaceProblem<dim>::solve()
  {
    LA::MPI::Vector completely_distributed_solution(locally_owned_dofs,
                                                    mpi_communicator);

    SolverControl solver_control(dof_handler.n_dofs(), 1e-12);

#ifdef USE_PETSC_LA
    // LA::SolverCG solver(solver_control, mpi_communicator);
    LA::SolverCG solver(solver_control);
#else
    LA::SolverCG solver(solver_control);
#endif

    LA::MPI::PreconditionAMG                 preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
    data.symmetric_operator = true;
#endif

    preconditioner.initialize(system_matrix, data);

    check_solver_within_range(solver.solve(system_matrix,
                                           completely_distributed_solution,
                                           system_rhs,
                                           preconditioner),
                              solver_control.last_step(),
                              1,
                              9);

    constraints.distribute(completely_distributed_solution);
    locally_relevant_solution = completely_distributed_solution;
  }

  template <int dim>
  void
  LaplaceProblem<dim>::output_results() const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution, "u");

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches(mapping);
    data_out.write_vtu_with_pvtu_record(
      "./", "solution", 0, mpi_communicator, 2, 8);
  }

  template <int dim>
  void
  LaplaceProblem<dim>::run()
  {
    deallog << "Running on "
            << Utilities::MPI::n_mpi_processes(mpi_communicator)
            << " MPI rank(s)..." << std::endl;
    GridIn<2, 2> grid_in;
    grid_in.attach_triangulation(triangulation);
    grid_in.read_partitioned_msh(SOURCE_DIR "/../grid/grids/unit-square");


    deallog << "Triangulation has " << triangulation.n_active_cells()
            << " active cells." << std::endl;

    deallog << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
            << ": Local cells = " << triangulation.n_active_cells()
            << ", Global cells = " << triangulation.n_global_active_cells()
            << std::endl;

    auto ghost_owners = triangulation.ghost_owners();
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": Ghost owners: ";
    for (const auto &owner : ghost_owners)
      std::cout << owner << " ";
    std::cout << std::endl;

    unsigned int n_locally_owned = 0;
    unsigned int n_ghost         = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          ++n_locally_owned;
        else if (cell->is_ghost())
          ++n_ghost;
      }

    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": Locally owned cells = " << n_locally_owned
              << ", Ghost cells = " << n_ghost << std::endl;

    MPI_Barrier(mpi_communicator);

    setup_system();

    deallog << "   Number of active cells:       "
            << triangulation.n_global_active_cells() << std::endl
            << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

    assemble_system();

    solve();

    if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
      output_results();
  }

} // namespace Step40

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // mpi_initlog();
  MPILogInitAll mpilog;


  try
    {
      using namespace Step40;
      LaplaceProblem<2> laplace_problem_2d;
      laplace_problem_2d.run();
    }
  catch (const std::exception &exc)
    {
      std::cerr << "\n\n----------------------------------------------------\n"
                << "Exception on processing: " << exc.what() << '\n'
                << "Aborting!\n"
                << "----------------------------------------------------\n"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << "\n\n----------------------------------------------------\n"
                << "Unknown exception!\n"
                << "Aborting!\n"
                << "----------------------------------------------------\n"
                << std::endl;
      return 1;
    }

  return 0;
}
