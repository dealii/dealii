/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2013 - 2021 by the deal.II authors
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
 * Author: Wolfgang Bangerth, Colorado State University, 2023
 */


#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/petsc_ts.h>

#include <fstream>
#include <iostream>


namespace Step86
{
  using namespace dealii;


  template <int dim>
  class HeatEquation : public ParameterAcceptor
  {
  public:
    HeatEquation(const MPI_Comm mpi_communicator = MPI_COMM_WORLD);
    void run();

  private:
    void setup_system();

    void implicit_function(const double                      time,
                           const PETScWrappers::MPI::Vector &solution,
                           const PETScWrappers::MPI::Vector &solution_dot,
                           PETScWrappers::MPI::Vector &      dst) const;

    void
    assemble_implicit_jacobian(const double                      time,
                               const PETScWrappers::MPI::Vector &solution,
                               const PETScWrappers::MPI::Vector &solution_dot,
                               const double                      shift);

    void distribute(const double time, PETScWrappers::MPI::Vector &dst) const;

    IndexSet algebraic_components() const;

    void output_results(const double                      time,
                        const PETScWrappers::MPI::Vector &solution,
                        const unsigned int timestep_number) const;

    void solve_with_jacobian(const PETScWrappers::MPI::Vector &src,
                             PETScWrappers::MPI::Vector &      dst) const;

    void
    prepare_for_coarsening_and_refinement(const PETScWrappers::MPI::Vector &y);

    void interpolate(const std::vector<PETScWrappers::MPI::Vector> &all_in,
                     std::vector<PETScWrappers::MPI::Vector> &      all_out);

    void fix_constraints(const double time) const;

    const MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim>                                 fe;
    DoFHandler<dim>                           dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double>         homogeneous_constraints;
    mutable AffineConstraints<double> constraints;

    PETScWrappers::MPI::SparseMatrix jacobian_matrix;

    PETScWrappers::MPI::Vector         solution;
    mutable PETScWrappers::MPI::Vector locally_relevant_solution;
    mutable PETScWrappers::MPI::Vector locally_relevant_solution_dot;

    PETScWrappers::TimeStepperData time_stepper_data;

    ParameterAcceptorProxy<Functions::ParsedFunction<dim>>
      initial_value_function;

    mutable ParameterAcceptorProxy<Functions::ParsedFunction<dim>>
      right_hand_side_function;

    mutable ParameterAcceptorProxy<Functions::ParsedFunction<dim>>
      boundary_values_function;

    unsigned int initial_global_refinement;
    unsigned int max_delta_refinement_level;
    unsigned int adaption_frequency;

    ConditionalOStream  pcout;
    mutable TimerOutput computing_timer;
  };


  template <int dim>
  HeatEquation<dim>::HeatEquation(const MPI_Comm mpi_communicator)
    : ParameterAcceptor("/Heat Equation/")
    , mpi_communicator(mpi_communicator)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe(1)
    , dof_handler(triangulation)
    , time_stepper_data("", "beuler", 0.0, 1.0, 0.025)
    , initial_value_function("/Heat Equation/Initial value")
    , right_hand_side_function("/Heat Equation/Right hand side")
    , boundary_values_function("/Heat Equation/Boundary values")
    , initial_global_refinement(5)
    , max_delta_refinement_level(2)
    , adaption_frequency(0)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
  {
    enter_subsection("Time stepper");
    enter_my_subsection(this->prm);
    time_stepper_data.add_parameters(this->prm);
    leave_my_subsection(this->prm);
    leave_subsection();

    add_parameter("Initial global refinement",
                  initial_global_refinement,
                  "Number of times the mesh is refined globally before "
                  "starting the time stepping.");
    add_parameter("Maximum delta refinement level",
                  max_delta_refinement_level,
                  "Maximum number of local refinement levels.");
    add_parameter("Adaption frequency",
                  adaption_frequency,
                  "When to adapt the mesh.");
  }



  template <int dim>
  void HeatEquation<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup_system");
    dof_handler.distribute_dofs(fe);

    pcout << std::endl
          << "===========================================" << std::endl
          << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl
          << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl
          << std::endl;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    homogeneous_constraints.clear();
    homogeneous_constraints.reinit(locally_relevant_dofs);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             homogeneous_constraints);
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            homogeneous_constraints);
    homogeneous_constraints.make_consistent_in_parallel(locally_owned_dofs,
                                                        locally_relevant_dofs,
                                                        mpi_communicator);
    homogeneous_constraints.close();

    fix_constraints(boundary_values_function.get_time());

    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    homogeneous_constraints,
                                    false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);

    // directly initialize from dsp, no need for the regular sparsity pattern:
    jacobian_matrix.reinit(locally_owned_dofs,
                           locally_owned_dofs,
                           dsp,
                           mpi_communicator);

    solution.reinit(locally_owned_dofs, mpi_communicator);
    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    locally_relevant_solution_dot.reinit(locally_owned_dofs,
                                         locally_relevant_dofs,
                                         mpi_communicator);
  }


  template <int dim>
  void HeatEquation<dim>::fix_constraints(const double time) const
  {
    TimerOutput::Scope t(computing_timer, "fix_constraints");
    boundary_values_function.set_time(time);
    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             boundary_values_function,
                                             constraints);
    constraints.make_consistent_in_parallel(locally_owned_dofs,
                                            locally_relevant_dofs,
                                            mpi_communicator);
    constraints.close();
  }


  template <int dim>
  void
  HeatEquation<dim>::implicit_function(const double                      time,
                                       const PETScWrappers::MPI::Vector &y,
                                       const PETScWrappers::MPI::Vector &y_dot,
                                       PETScWrappers::MPI::Vector &dst) const
  {
    TimerOutput::Scope t(computing_timer, "implicit_function");
    right_hand_side_function.set_time(time);

    PETScWrappers::MPI::Vector tmp_solution(y);
    PETScWrappers::MPI::Vector tmp_solution_dot(y_dot);

    // Fix boundary conditions
    fix_constraints(time);
    constraints.distribute(tmp_solution);
    homogeneous_constraints.distribute(tmp_solution_dot);

    // Copy to ghosted vectors
    locally_relevant_solution     = tmp_solution;
    locally_relevant_solution_dot = tmp_solution_dot;

    QGauss<dim>   quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> solution_gradients(n_q_points);
    std::vector<double>         solution_dot_values(n_q_points);

    Vector<double> cell_residual(dofs_per_cell);

    dst = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          fe_values.get_function_gradients(locally_relevant_solution,
                                           solution_gradients);
          fe_values.get_function_values(locally_relevant_solution_dot,
                                        solution_dot_values);

          cell->get_dof_indices(local_dof_indices);

          cell_residual = 0;
          for (const unsigned int q : fe_values.quadrature_point_indices())
            for (const unsigned int i : fe_values.dof_indices())
              {
                cell_residual(i) +=
                  (fe_values.shape_value(i, q) * solution_dot_values[q] +
                   fe_values.shape_grad(i, q) * solution_gradients[q] -
                   fe_values.shape_value(i, q) *
                     right_hand_side_function.value(
                       fe_values.quadrature_point(q))) *
                  fe_values.JxW(q);
              }
          homogeneous_constraints.distribute_local_to_global(cell_residual,
                                                             local_dof_indices,
                                                             dst);
        }
    dst.compress(VectorOperation::add);

    // Now we correct the entries corresponding to constrained degrees of
    // freedom. We only care about essential boundaries where tmp_solution
    // already contains the right values.
    for (const auto &c : homogeneous_constraints.get_lines())
      if (locally_owned_dofs.is_element(c.index))
        {
          if (c.entries.empty())
            {
              // This is an essential boundary condition
              dst[c.index] = y[c.index] - tmp_solution[c.index];
            }
          else
            {
              // This is a hanging node. We set it equal to y because we want
              // an exact Jacobian. We ignore hanging nodes in the computation
              // of the local truncation error, since we always call distribute
              // after a successful stage.
              dst[c.index] = y[c.index];
            }
        }
    dst.compress(VectorOperation::insert);
  }


  template <int dim>
  void HeatEquation<dim>::assemble_implicit_jacobian(
    const double,
    const PETScWrappers::MPI::Vector &,
    const PETScWrappers::MPI::Vector &,
    const double shift)
  {
    TimerOutput::Scope t(computing_timer, "assemble_implicit_jacobian");
    QGauss<dim>        quadrature_formula(fe.degree + 1);
    FEValues<dim>      fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    jacobian_matrix = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          cell->get_dof_indices(local_dof_indices);

          cell_matrix = 0;
          for (const unsigned int q : fe_values.quadrature_point_indices())
            for (const unsigned int i : fe_values.dof_indices())
              for (const unsigned int j : fe_values.dof_indices())
                {
                  cell_matrix(i, j) +=
                    (shift * fe_values.shape_value(i, q) *
                       fe_values.shape_value(j, q) +
                     fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)) *
                    fe_values.JxW(q);
                }
          homogeneous_constraints.distribute_local_to_global(cell_matrix,
                                                             local_dof_indices,
                                                             jacobian_matrix);
        }
    jacobian_matrix.compress(VectorOperation::add);
    // Now we correct the entries corresponding to constrained degrees of
    // freedom. We want the Jacobian to be one on constrained dofs
    for (const auto &c : homogeneous_constraints.get_lines())
      jacobian_matrix.set(c.index, c.index, 1.0);
    jacobian_matrix.compress(VectorOperation::insert);
  }

  template <int dim>
  void HeatEquation<dim>::distribute(const double                time,
                                     PETScWrappers::MPI::Vector &dst) const
  {
    TimerOutput::Scope t(computing_timer, "distribute");
    fix_constraints(time);
    constraints.distribute(dst);
  }

  template <int dim>
  void
  HeatEquation<dim>::solve_with_jacobian(const PETScWrappers::MPI::Vector &src,
                                         PETScWrappers::MPI::Vector &dst) const
  {
    TimerOutput::Scope      t(computing_timer, "solve_with_jacobian");
    SolverControl           solver_control(1000, 1e-8 * src.l2_norm());
    PETScWrappers::SolverCG cg(solver_control);
    cg.set_prefix("user_");

#if defined(PETSC_HAVE_HYPRE)
    PETScWrappers::PreconditionBoomerAMG preconditioner;
    preconditioner.initialize(jacobian_matrix);
#else
    PETScWrappers::PreconditionSSOR preconditioner;
    preconditioner.initialize(
      jacobian_matrix, PETScWrappers::PreconditionSSOR::AdditionalData(1.0));
#endif
    cg.solve(jacobian_matrix, dst, src, preconditioner);

    pcout << "     " << solver_control.last_step() << " linear iterations."
          << std::endl;
  }



  template <int dim>
  void
  HeatEquation<dim>::output_results(const double                      time,
                                    const PETScWrappers::MPI::Vector &y,
                                    const unsigned int timestep_number) const
  {
    TimerOutput::Scope t(computing_timer, "output_results");
    DataOut<dim>       data_out;

    // PETScWrappers::MPI::Vector tmp_solution(y);
    // fix_constraints(time);
    // constraints.distribute(tmp_solution);
    locally_relevant_solution = y;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution, "U");

    data_out.build_patches();

    data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number));

    const std::string filename =
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu";
    data_out.write_vtu_in_parallel(filename, mpi_communicator);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        static std::vector<std::pair<double, std::string>> times_and_names;
        times_and_names.emplace_back(time, filename);
        std::ofstream pvd_output("solution.pvd");
        DataOutBase::write_pvd_record(pvd_output, times_and_names);
      }
  }


  template <int dim>
  void HeatEquation<dim>::prepare_for_coarsening_and_refinement(
    const PETScWrappers::MPI::Vector &sol)
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    locally_relevant_solution = sol;
    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1),
                                       {},
                                       locally_relevant_solution,
                                       estimated_error_per_cell);

    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
      triangulation, estimated_error_per_cell, 0.6, 0.4);

    const auto max_grid_level =
      initial_global_refinement + max_delta_refinement_level;
    const auto min_grid_level = initial_global_refinement;

    if (triangulation.n_levels() > max_grid_level)
      for (const auto &cell :
           triangulation.active_cell_iterators_on_level(max_grid_level))
        cell->clear_refine_flag();
    for (const auto &cell :
         triangulation.active_cell_iterators_on_level(min_grid_level))
      cell->clear_coarsen_flag();
  }



  template <int dim>
  void HeatEquation<dim>::interpolate(
    const std::vector<PETScWrappers::MPI::Vector> &all_in,
    std::vector<PETScWrappers::MPI::Vector> &      all_out)
  {
    parallel::distributed::SolutionTransfer<dim, PETScWrappers::MPI::Vector>
      solution_trans(dof_handler);

    std::vector<PETScWrappers::MPI::Vector> all_in_ghosted(all_in.size());
    std::vector<const PETScWrappers::MPI::Vector *> all_in_ghosted_ptr(
      all_in.size());
    std::vector<PETScWrappers::MPI::Vector *> all_out_ptr(all_in.size());
    for (unsigned int i = 0; i < all_in.size(); ++i)
      {
        all_in_ghosted[i].reinit(locally_owned_dofs,
                                 locally_relevant_dofs,
                                 mpi_communicator);
        all_in_ghosted[i]     = all_in[i];
        all_in_ghosted_ptr[i] = &all_in_ghosted[i];
      }

    triangulation.prepare_coarsening_and_refinement();
    solution_trans.prepare_for_coarsening_and_refinement(all_in_ghosted_ptr);
    triangulation.execute_coarsening_and_refinement();

    setup_system();

    all_out.resize(all_in.size());
    for (unsigned int i = 0; i < all_in.size(); ++i)
      {
        all_out[i].reinit(locally_owned_dofs, mpi_communicator);
        all_out_ptr[i] = &all_out[i];
      }
    solution_trans.interpolate(all_out_ptr);
  }



  template <int dim>
  IndexSet HeatEquation<dim>::algebraic_components() const
  {
    return DoFTools::extract_hanging_node_dofs(dof_handler);
  }



  template <int dim>
  void HeatEquation<dim>::run()
  {
    GridGenerator::hyper_L(triangulation);
    triangulation.refine_global(initial_global_refinement);

    setup_system();

    VectorTools::interpolate(dof_handler, initial_value_function, solution);

    PETScWrappers::TimeStepper<PETScWrappers::MPI::Vector,
                               PETScWrappers::MPI::SparseMatrix>
      petsc_ts(time_stepper_data);

    petsc_ts.set_matrices(jacobian_matrix, jacobian_matrix);

    petsc_ts.implicit_function =
      [&](const auto t, const auto &y, const auto &y_dot, auto &res) {
        this->implicit_function(t, y, y_dot, res);
      };

    petsc_ts.setup_jacobian =
      [&](const auto t, const auto &y, const auto &y_dot, const auto alpha) {
        this->assemble_implicit_jacobian(t, y, y_dot, alpha);
      };

    petsc_ts.solve_with_jacobian = [&](const auto &src, auto &dst) {
      this->solve_with_jacobian(src, dst);
    };

    petsc_ts.distribute = [&](const auto t, auto &y) {
      this->distribute(t, y);
    };

    petsc_ts.algebraic_components = [&]() {
      return this->algebraic_components();
    };

    petsc_ts.monitor =
      [&](const auto t, const auto &y, const auto step_number) {
        pcout << "Time step " << step_number << " at t=" << t << std::endl;
        this->output_results(t, y, step_number);
      };

    petsc_ts.prepare_for_coarsening_and_refinement =
      [&](const auto, const auto it, const auto &y, auto &resize) {
        if (it > 0 && this->adaption_frequency > 0 &&
            it % this->adaption_frequency == 0)
          {
            this->prepare_for_coarsening_and_refinement(y);
            resize = true;
          }
      };

    petsc_ts.interpolate = [&](const auto &all_in, auto &all_out) {
      this->interpolate(all_in, all_out);
      petsc_ts.set_matrices(this->jacobian_matrix, this->jacobian_matrix);
    };

    petsc_ts.solve(solution);
  }
} // namespace Step86



int main(int argc, char **argv)
{
  try
    {
      using namespace Step86;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      HeatEquation<2>                  heat_equation_solver;
      ParameterAcceptor::initialize("heat_equation.prm",
                                    "heat_equation_used.prm");
      heat_equation_solver.run();
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
