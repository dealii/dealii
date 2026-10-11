/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Author: David Wells, University of North Carolina, Chapel Hill, 2021
 */

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools_boundary.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_rhs.h>

#include <deal.II/sundials/arkode.h>

#include <iostream>
#include <memory>

#define ODE 0
constexpr int fe_degree = 1;

using namespace dealii;

class ExactSolution : public Function<2>
{
  public:
    virtual double value(const Point<2> &point, const unsigned int /*component*/ = 0) const override
    {
      // return std::sin(point[0])*std::sin(point[1])* get_time();
      const double &x = point[0];
      const double &y = point[1];
      return x * (1 - x) * y * (1 - y) * get_time();
    }
};

class Forcing : public Function<2>
{
  public:
    virtual double value(const Point<2> &point, const unsigned int /*component*/ = 0) const override
    {
      const double &x = point[0];
      const double &y = point[1];
#if ODE
      // return std::sin(point[0]) * std::sin(point[1]);
      return x * (1 - x) * y * (1 - y);
#else
      // return (2 * get_time() + 1.0) * std::sin(point[0]) * std::sin(point[1]);
      return x * (x - 1) * y * (y - 1) - 2 * get_time() * (x * (x - 1) + y * (y - 1));
#endif
    }
};

namespace Step84
{
  class QGESolver
  {
    public:
      QGESolver(const unsigned int n_refinements);

      double run();

      // Today's goal - heat equation

    private:
      ConditionalOStream pcout;
      TimerOutput        computing_timer;

      unsigned int n_refinements;

      Triangulation<2> triangulation;
      FE_Q<2> fe;
      QGauss<2> quadrature;
      QGauss<2> rhs_quadrature;
      DoFHandler<2> dof_handler;

      std::shared_ptr<Utilities::MPI::Partitioner> partitioner;

      AffineConstraints<double> homogeneous_constraints;

      MappingCartesian<2> mapping;
  };

  QGESolver::QGESolver(const unsigned int n_refinements)
    : pcout(std::cout,
            (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
    , computing_timer(MPI_COMM_WORLD,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , n_refinements(n_refinements)
    , fe(fe_degree)
    , quadrature(fe_degree + 1)
    , rhs_quadrature(fe_degree + 3)
    , dof_handler(triangulation)
  {
    TimerOutput::Scope t(computing_timer, "setup grid and dofs");
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(n_refinements);
    dof_handler.distribute_dofs(fe);

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    partitioner = std::make_shared<Utilities::MPI::Partitioner>(
      dof_handler.locally_owned_dofs(),
      locally_relevant_dofs,
      triangulation.get_communicator());

    // This only works with time-independent constraints since we don't have
    // access to the current time in mass solver
    VectorTools::interpolate_boundary_values(mapping,
                                             dof_handler,
                                             0,
                                             ZeroFunction<2>(),
                                             homogeneous_constraints);
    homogeneous_constraints.close();

    pcout << "Number of dofs = " << dof_handler.n_dofs() << std::endl;
  }

  double
  QGESolver::run()
  {
    TimerOutput::Scope t1(computing_timer, "setup operators");

    using VectorType = LinearAlgebra::distributed::Vector<double>;
    SUNDIALS::ARKode<VectorType>::AdditionalData arkode_parameters;
    // do one time step:
    arkode_parameters.final_time = 1.0;
    arkode_parameters.initial_step_size = triangulation.begin_active()->diameter() / 4.0;

    // it's the heat equation:
    arkode_parameters.implicit_function_is_linear = true;
    arkode_parameters.implicit_function_is_time_independent = true;
    arkode_parameters.mass_is_time_independent = true;

    arkode_parameters.maximum_order = fe_degree + 1;

    // Set up linear operators:
    auto matrix_free = std::make_shared<MatrixFree<2, double>>();
    matrix_free->reinit(mapping, dof_handler, homogeneous_constraints, quadrature);

    MatrixFreeOperators::MassOperator<2, fe_degree, fe_degree + 1, 1> mass;
    mass.initialize(matrix_free);
    mass.compute_diagonal();
    PreconditionJacobi<MatrixFreeOperators::Base<2>> mass_preconditioner;
    mass_preconditioner.initialize(mass, 1.0);

    MatrixFreeOperators::LaplaceOperator<2, fe_degree, fe_degree + 1, 1> laplace;
    laplace.initialize(matrix_free);
    laplace.compute_diagonal();
    PreconditionJacobi<MatrixFreeOperators::Base<2>> laplace_preconditioner;
    laplace_preconditioner.initialize(laplace, 1.0);
    t1.stop();

    SUNDIALS::ARKode<VectorType>
      time_integrator(arkode_parameters, triangulation.get_communicator());

#if ODE
    // no implicitly integrated function
#else
    time_integrator.implicit_function = [&](const double /*time*/,
                                            const VectorType &solution,
                                            VectorType &implicit_result)
    {
      TimerOutput::Scope t1(computing_timer, "implicit function");
      // TODO - the built-in Laplace operator is positive definite
      laplace.vmult(implicit_result, solution);
      implicit_result *= -1.0;

      return 0;
    };
#endif

    // Forcing is explicit - it doesn't create a time-step restriction
    time_integrator.explicit_function = [&](const double time,
                                            const VectorType &/*solution*/,
                                            VectorType &explicit_result)
    {
      TimerOutput::Scope t1(computing_timer, "explicit function");
      explicit_result = 0.0;

      Forcing forcing;
      forcing.set_time(time);
      VectorTools::create_right_hand_side(mapping,
                                          dof_handler,
                                          rhs_quadrature,
                                          forcing,
                                          explicit_result,
                                          homogeneous_constraints);

      return 0;
    };

    time_integrator.mass_times_vector = [&](const double /*time*/,
                                            const VectorType &solution,
                                            VectorType &mass_result)
    {
      TimerOutput::Scope t1(computing_timer, "mass times vector");
      // NOTE - this isn't enough since we have constraints. Note that somewhere.
      mass.vmult(mass_result, solution);

      return 0;
    };

    time_integrator.jacobian_times_vector = [&](const VectorType &src,
                                                VectorType &dst,
                                                const double time,
                                                const VectorType &/*solution*/,
                                                const VectorType &/*implicit_rhs*/)
    {
      // TODO - the wrappers aren't smart enough to figure out that this should
      // be done if the implicit function is linear.
      return time_integrator.implicit_function(time, src, dst);
    };


    // TODO - we need the current time to enforce time-dependent constraints
    // correctly. We either need a backdoor (there may be one - I didn't look
    // very hard) or we should encode it in SundialsOperator.
    time_integrator.solve_mass = [&](SUNDIALS::SundialsOperator<VectorType> &op,
                                     SUNDIALS::SundialsPreconditioner<VectorType> &prec,
                                     VectorType &x,
                                     const VectorType &b,
                                     double /*tolerance*/)
    {
      TimerOutput::Scope t1(computing_timer, "solve mass");
      // TODO - is tolerance relative or absolute?
      SolverControl control(1000, 1e-10 * b.l2_norm()); // TODO - use real tolerance
      SolverCG<VectorType> cg(control);

      cg.solve(op, x, b, prec);
      homogeneous_constraints.distribute(x);

      return 0;
    };

    // TODO - what we really want to use here is a proper multigrid
    // preconditioner for a Helmholtz-like problem M - gamma J, but we don't
    // know gamma so we cannot set this up ourselves. This is a limitation in
    // the wrapper (the pre-4.0 versions passed gamma explicitly).
    //
    // We can use ARKStepGetCurrentGamma() if necessary
    time_integrator.solve_linearized_system = [&](
      SUNDIALS::SundialsOperator<VectorType> &op,
      SUNDIALS::SundialsPreconditioner<VectorType> &prec,
      VectorType &x,
      const VectorType &b,
      double /*tolerance*/)
    {
      TimerOutput::Scope t1(computing_timer, "solve linearized");
      // TODO - same notes on tolerance
      SolverControl control(1000, 1e-10 * b.l2_norm());
      SolverGMRES<VectorType> solver(control);

      solver.solve(op, x, b, prec);
      homogeneous_constraints.distribute(x);

      return 0;
    };

    VectorType solution(partitioner);
    const double time_step = triangulation.begin_active()->diameter()/2.0;
    unsigned int step_n = 0;
    double next_time = 0.0;
    do
    {
      ++step_n;
      next_time = std::min(step_n * time_step, arkode_parameters.final_time);
      time_integrator.solve_ode_incrementally(solution, next_time);
    } while (next_time < arkode_parameters.final_time);

    ExactSolution exact_solution;
    exact_solution.set_time(arkode_parameters.final_time);
    Vector<double> cell_errors(triangulation.n_active_cells());
    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      exact_solution,
                                      cell_errors,
                                      rhs_quadrature,
                                      VectorTools::L2_norm);

    const double global_error =
    VectorTools::compute_global_error(triangulation,
                                      cell_errors,
                                      VectorTools::L2_norm);
    pcout << "Number of macro steps = " << step_n << std::endl;
    pcout << "Global error = " << global_error << std::endl;

    {
      DataOut<2> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "u");
      data_out.add_data_vector(cell_errors, "e");
      data_out.build_patches();

      data_out.write_vtu_with_pvtu_record(
        "./", "solution", n_refinements, triangulation.get_communicator(), 2, 8);

    }

    return global_error;
  }
}

// @sect3{The <code>main</code> Function}

// Having made it this far, there is, again, nothing much to discuss for the
// main function of this program: it looks like all such functions since step-6.
int main(int argc, char **argv)
{
  try
    {
      using namespace Step84;
      using namespace dealii;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      ConditionalOStream pcout(std::cout,
                               (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

      double old_error = 1.0;
      for (unsigned int n_refinements = 2; n_refinements < 8; ++n_refinements)
        {
          QGESolver solver(n_refinements);
          double new_error = solver.run();

          pcout << "ratio = " << old_error / new_error << std::endl;
          old_error = new_error;
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
