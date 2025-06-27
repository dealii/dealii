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



// solves a 2D Poisson equation for linear elements with the simple
// (Ifpack-based) Trilinos preconditioners

#include "deal.II/base/exceptions.h"
#include "deal.II/base/logstream.h"
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "deal.II/lac/trilinos_tpetra_vector.h"
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
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

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dof_handler;

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
  value(const Point<dim> &p, const unsigned int component = 0) const;
};



template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;
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
  : fe(1)
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

  constraints.clear();
  std::map<unsigned int, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           BoundaryValues<dim>(),
                                           constraints);
  constraints.close();

  DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, c_sparsity, constraints, false);
  system_matrix.reinit(c_sparsity);
  IndexSet is(dof_handler.n_dofs());
  is.add_range(0, dof_handler.n_dofs());
  solution.reinit(is);
  system_rhs.reinit(is);
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

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();

  for (; cell != endc; ++cell)
    {
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
    deallog.push("Identity");
    static constexpr std::array<unsigned int, 2> lower{{49, 100}};
    LinearAlgebra::TpetraWrappers::PreconditionIdentity<double,
                                                        MemorySpace::Default>
      preconditioner;
    solution = 0;
    SolverControl solver_control(1000, 1e-10);
    SolverCG<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
                           solver(solver_control);
    Teuchos::ParameterList precon_params;
    preconditioner.initialize(system_matrix);
    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      lower[cycle],
      lower[cycle] + 2);
    deallog.pop();
  }

  {
    // Make sure the general Ifpack preconditioner works at all.
    // For this we use the point relaxation (SSOR)
    // with default parameters, resulting in symmetric Gauss-Seidel (SGS)
    deallog.push("CustomSGS");
    static constexpr std::array<unsigned int, 2> lower{{49, 100}};
    LinearAlgebra::TpetraWrappers::PreconditionIfpack<double,
                                                      MemorySpace::Default>
      preconditioner("RELAXATION");
    solution = 0;
    SolverControl solver_control(1000, 1e-10);
    SolverCG<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
                           solver(solver_control);
    Teuchos::ParameterList precon_params;
    precon_params.set("relaxation: type", "Symmetric Gauss-Seidel");
    preconditioner.set_parameter_list(precon_params);
    preconditioner.initialize(system_matrix);
    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      lower[cycle],
      lower[cycle] + 2);
    deallog.pop();
  }

  {
    deallog.push("Jacobi");
    static constexpr std::array<unsigned int, 2> lower{{49, 100}};
    LinearAlgebra::TpetraWrappers::PreconditionJacobi<double,
                                                      MemorySpace::Default>
      preconditioner;
    solution = 0;

    SolverControl solver_control(1000, 1e-10);
    SolverCG<
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
    deallog.push("l1Jacobi");
    static constexpr std::array<unsigned int, 2> lower{{49, 100}};
    LinearAlgebra::TpetraWrappers::PreconditionL1Jacobi<double,
                                                        MemorySpace::Default>
      preconditioner;
    solution = 0;

    SolverControl solver_control(1000, 1e-10);
    SolverCG<
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
    deallog.push("l1GaussSeidel");
    static constexpr std::array<unsigned int, 2> lower{{31, 62}};
    LinearAlgebra::TpetraWrappers::
      PreconditionL1GaussSeidel<double, MemorySpace::Default>
        preconditioner;
    solution = 0;

    SolverControl solver_control(1000, 1e-10);
    SolverBicgstab<
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
    deallog.push("SOR");
    static constexpr std::array<unsigned int, 2> lower{{31, 62}};
    LinearAlgebra::TpetraWrappers::PreconditionSOR<double, MemorySpace::Default>
      preconditioner;
    solution = 0;

    SolverControl solver_control(1000, 1e-5);
    SolverBicgstab<
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
    deallog.push("SSOR");
    static constexpr std::array<unsigned int, 2> lower{{40, 77}};
    LinearAlgebra::TpetraWrappers::PreconditionSSOR<double,
                                                    MemorySpace::Default>
      preconditioner;
    solution = 0;

    SolverControl solver_control(1000, 1e-10);
    SolverCG<
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
    deallog.push("BlockJacobi");
    static constexpr std::array<unsigned int, 2> lower{{73, 145}};
    LinearAlgebra::TpetraWrappers::PreconditionBlockJacobi<double,
                                                           MemorySpace::Default>
      preconditioner;
    LinearAlgebra::TpetraWrappers::PreconditionBlockJacobi<
      double,
      MemorySpace::Default>::AdditionalData data;
    data.n_local_parts = 16;
    solution           = 0;
    SolverControl solver_control(1000, 1e-10);
    SolverCG<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
      solver(solver_control);
    preconditioner.initialize(system_matrix, data);
    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      lower[cycle],
      lower[cycle] + 2);
    deallog.pop();
  }

  {
    deallog.push("BlockSOR");
    static constexpr std::array<unsigned int, 2> lower{{18, 37}};
    LinearAlgebra::TpetraWrappers::PreconditionBlockSOR<double,
                                                        MemorySpace::Default>
      preconditioner;
    LinearAlgebra::TpetraWrappers::
      PreconditionBlockSOR<double, MemorySpace::Default>::AdditionalData data;
    data.n_local_parts = 16;
    data.omega         = 0.8;
    solution           = 0;
    SolverControl solver_control(1000, 1e-5);
    SolverBicgstab<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
      solver(solver_control);
    preconditioner.initialize(system_matrix, data);
    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      lower[cycle],
      lower[cycle] + 2);
    deallog.pop();
  }

  {
    deallog.push("BlockSSOR");
    static constexpr std::array<unsigned int, 2> lower{{30, 59}};
    LinearAlgebra::TpetraWrappers::PreconditionBlockSSOR<double,
                                                         MemorySpace::Default>
      preconditioner;
    LinearAlgebra::TpetraWrappers::
      PreconditionBlockSSOR<double, MemorySpace::Default>::AdditionalData data;
    data.n_local_parts = 16;
    data.omega         = 1.2;
    solution           = 0;
    SolverControl solver_control(1000, 1e-10);
    SolverCG<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
      solver(solver_control);
    preconditioner.initialize(system_matrix, data);
    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      lower[cycle],
      lower[cycle] + 2);
    deallog.pop();
  }

  {
    static constexpr std::array<unsigned int, 2> lower{{30, 56}};
    deallog.push("ILU");
    LinearAlgebra::TpetraWrappers::PreconditionILU<double, MemorySpace::Default>
      preconditioner;
    solution = 0;
    SolverControl solver_control(1000, 1e-10);
    SolverCG<
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
    deallog.push("ILUT");
    static constexpr std::array<unsigned int, 2> lower{{11, 19}};
    LinearAlgebra::TpetraWrappers::PreconditionILUT<double,
                                                    MemorySpace::Default>
      preconditioner;
    LinearAlgebra::TpetraWrappers::
      PreconditionILUT<double, MemorySpace::Default>::AdditionalData data;
    data.ilut_drop = 1e-6;
    data.ilut_fill = 3;
    solution       = 0;
    SolverControl solver_control(1000, 1e-5);
    SolverBicgstab<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
      solver(solver_control);
    preconditioner.initialize(system_matrix, data);
    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      lower[cycle],
      lower[cycle] + 2);
    deallog.pop();
  }

  {
    deallog.push("Chebyshev");
    static constexpr std::array<unsigned int, 2> lower{{23, 46}};
    LinearAlgebra::TpetraWrappers::PreconditionChebyshev<double,
                                                         MemorySpace::Default>
      preconditioner;
    LinearAlgebra::TpetraWrappers::
      PreconditionChebyshev<double, MemorySpace::Default>::AdditionalData data;
    data.max_eigenvalue = 2.5;
    data.degree         = 3;
    solution            = 0;
    SolverControl solver_control(1000, 1e-10);
    SolverCG<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
      solver(solver_control);
    preconditioner.initialize(system_matrix, data);
    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      lower[cycle],
      lower[cycle] + 2);
    deallog.pop();
  }

  //   {
  //     deallog.push("Direct");
  //     TrilinosWrappers::PreconditionBlockwiseDirect preconditioner;
  //     solution = 0;
  //     SolverControl solver_control(1000, 1e-10);
  //     SolverCG<>    solver(solver_control);
  //     preconditioner.initialize(system_matrix);
  //     check_solver_within_range(
  //       solver.solve(system_matrix, solution, system_rhs, preconditioner),
  //       solver_control.last_step(),
  //       1U,
  //       1U);
  //     deallog.pop();
  //   }

  {
    // Make sure the general Ifpack preconditioner throws an Exception for an
    // unsupported solver.
    // IGLOO (incomplete global least optimal option) is of course a made up
    // preconditioner name.

    deallog.push("CustomIGLOO");
    try
      {
        LinearAlgebra::TpetraWrappers::PreconditionIfpack<double,
                                                          MemorySpace::Default>
          preconditioner("IGLOO");
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
