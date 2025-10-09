// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// solves a 2D Poisson equation for linear elements with the FROSch
// preconditioners

#include "deal.II/base/exceptions.h"
#include "deal.II/base/logstream.h"
#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include "deal.II/lac/trilinos_tpetra_vector.h"
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_tpetra_precondition.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include <Teuchos_ParameterList.hpp>

#include "../tests.h"



template <int dim>
class Step4
{
public:
  Step4();
  void
  run();

private:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve(int cycle);

  MPI_Comm mpi_communicator;

  parallel::distributed::Triangulation<dim> triangulation;
  FE_Q<dim>                                 fe;
  DoFHandler<dim>                           dof_handler;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  AffineConstraints<double> constraints;

  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default>
    system_matrix;

  LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> solution;
  LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>
    system_rhs;
};


template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};



template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};



template <int dim>
double
RightHandSide<dim>::value(const Point<dim> &p,
                          const unsigned int /*component*/) const
{
  double return_value = 0;
  for (unsigned int i = 0; i < dim; ++i)
    return_value += 4 * std::pow(p[i], 4);

  return return_value;
}



template <int dim>
double
BoundaryValues<dim>::value(const Point<dim> &p,
                           const unsigned int /*component*/) const
{
  return p.square();
}



template <int dim>
Step4<dim>::Step4()
  : mpi_communicator(MPI_COMM_WORLD)
  , triangulation(mpi_communicator)
  , fe(1)
  , dof_handler(triangulation)
{}


template <int dim>
void
Step4<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(5);
}



template <int dim>
void
Step4<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  std::map<unsigned int, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           BoundaryValues<dim>(),
                                           constraints);
  constraints.close();

  // Note: This is a completely distributed solution vector.
  solution.reinit(locally_owned_dofs, mpi_communicator);

  system_rhs.reinit(locally_owned_dofs,
                    locally_relevant_dofs,
                    mpi_communicator,
                    true);

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
Step4<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  const RightHandSide<dim> right_hand_side;

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_point) *
                 fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));

            cell_rhs(i) +=
              (fe_values.shape_value(i, q_point) *
               right_hand_side.value(fe_values.quadrature_point(q_point)) *
               fe_values.JxW(q_point));
          }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
  system_matrix.compress(VectorOperation::add);
}



template <int dim>
void
Step4<dim>::solve(int cycle)
{
  deallog.push(Utilities::int_to_string(dof_handler.n_dofs(), 5));

  {
    deallog.push("FROSch - one level");
    // Note: The numer of iteration depends on the number of MPI-ranks.
    static constexpr std::array<unsigned int, 2> lower{{15, 23}};
    LinearAlgebra::TpetraWrappers::PreconditionFROSch<double,
                                                      MemorySpace::Default>
      preconditioner("One Level");
    solution = 0;
    SolverControl solver_control(1000, 1e-10);
    SolverGMRES<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
      solver(solver_control);
    preconditioner.initialize(system_matrix);
    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      lower[cycle],
      lower[cycle] + 2);
    deallog.pop();
  }

  {
    deallog.push("FROSch - two level");
    // Note: The numer of iteration depends on the number of MPI-ranks.
    static constexpr std::array<unsigned int, 2> lower{{7, 11}};
    LinearAlgebra::TpetraWrappers::PreconditionFROSch<double,
                                                      MemorySpace::Default>
      preconditioner("One Level");
    solution = 0;
    SolverControl solver_control(1000, 1e-10);
    SolverGMRES<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
      solver(solver_control);
    preconditioner.initialize(system_matrix, dof_handler);
    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      lower[cycle],
      lower[cycle] + 2);
    deallog.pop();
  }

  {
    // Make sure the FROSch preconditioner throws an Exception for an
    // unsupported precondition method.
    // RIBBIT ('Frosch' is the german word for 'frog', and 'ribbit' is the
    // sound the frog makes) is of course a made up precondition method name.
    deallog.push("FROSch - RIBBIT");
    try
      {
        LinearAlgebra::TpetraWrappers::PreconditionFROSch<double,
                                                          MemorySpace::Default>
          preconditioner("RIBBIT");
      }
    catch (dealii::ExceptionBase &exc)
      {
        deallog << "Error: " << std::endl;
        exc.print_info(deallog.get_file_stream());
      }
    deallog.pop();
  }

  deallog.pop();
}



template <int dim>
void
Step4<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 2; ++cycle)
    {
      if (cycle == 0)
        make_grid();
      else
        triangulation.refine_global(1);

      setup_system();
      assemble_system();
      solve(cycle);
    }
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  try
    {
      Step4<2> test;
      test.run();
    }
  catch (const std::exception &exc)
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
