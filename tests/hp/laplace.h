// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// base header for hp-FEM test on Laplace equation.


#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_series.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/smoothness_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>

#include "../tests.h"



/**
 * Basic class for Laplace problem
 */
template <int dim>
class Laplace
{
public:
  Laplace(const Function<dim> &force_function,
          const Function<dim> &exact_solution,
          const Function<dim> &boundary_conditions,
          const unsigned int   n_cycles,
          const std::string    txt_file_name);

  virtual ~Laplace();

  void
  run();

  DoFHandler<dim> &
  get_dof_handler();

  void
  setup_solve_estimate(Vector<float> &output_estimate);

protected:
  void
  setup_system();

  virtual void
  setup_geometry() = 0;

  void
  assemble();

  virtual void
  solve();

  /**
   * estimate error
   */
  virtual void
  estimate_error() = 0;

  /**
   * mark cells for h-refinement based on error estimation only
   */
  virtual void
  mark_h_cells() = 0;

  /**
   * remove h-refinement flag from some cells and flag cells for p-refinement
   */
  virtual void
  substitute_h_for_p() = 0;

  void
  refine_grid(const unsigned int cycle);

  void
  calculate_error();

  void
  output_results(int cycle);

  void
  print_errors();

  const Function<dim> &force_function;
  const Function<dim> &exact_solution;
  const Function<dim> &boundary_conditions;

  Triangulation<dim>    triangulation;
  hp::FECollection<dim> fe;
  DoFHandler<dim>       dof_handler;
  hp::QCollection<dim>  quadrature;
  hp::QCollection<dim>  quadrature_infty;

  AffineConstraints<double>      constraints;
  SparsityPattern                sparsity_pattern;
  TrilinosWrappers::SparseMatrix system_matrix;
  TrilinosWrappers::MPI::Vector  system_rhs;
  TrilinosWrappers::MPI::Vector  solution;
  TrilinosWrappers::MPI::Vector  solution_locally_relevant;

  Vector<float> estimated_error_per_cell;
  double        total_error;

  std::pair<unsigned int, unsigned int> hp_number;

  double L2_error;
  double H1_error;
  double Linfty_error;

  const unsigned int n_cycles;

  std::string      sp;
  ConvergenceTable error_table;
  std::string      output_name;

  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;
  ConditionalOStream pcout;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
};



// implementation
template <int dim>
Laplace<dim>::Laplace(const Function<dim> &force_function,
                      const Function<dim> &exact_solution,
                      const Function<dim> &boundary_conditions,
                      const unsigned int   n_cycles,
                      const std::string    output_name)
  : force_function(force_function)
  , exact_solution(exact_solution)
  , boundary_conditions(boundary_conditions)
  , dof_handler(triangulation)
  , n_cycles(n_cycles)
  , sp(" ")
  , output_name(output_name)
  , mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
{
  hp_number.first  = 0;
  hp_number.second = 0;

  deallog << std::endl;
}



template <int dim>
Laplace<dim>::~Laplace()
{
  dof_handler.clear();
}



template <int dim>
DoFHandler<dim> &
Laplace<dim>::get_dof_handler()
{
  return dof_handler;
}



template <int dim>
void
Laplace<dim>::setup_solve_estimate(Vector<float> &output_estimate)
{
  setup_system();
  assemble();
  solve();
  estimate_error();
  output_estimate = estimated_error_per_cell;
}



template <int dim>
void
Laplace<dim>::setup_system()
{
  GridTools::partition_triangulation(n_mpi_processes, triangulation);

  dof_handler.distribute_dofs(fe);

  locally_owned_dofs = dof_handler.locally_owned_dofs();

  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  AssertThrow(locally_relevant_dofs.n_elements() == dof_handler.n_dofs(),
              ExcInternalError());

  // init vectors
  solution_locally_relevant.reinit(locally_owned_dofs,
                                   locally_relevant_dofs,
                                   mpi_communicator);
  solution.reinit(locally_owned_dofs, mpi_communicator);
  solution                  = 0;
  solution_locally_relevant = solution;
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs = 0;

  constraints.clear();
  // constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           boundary_conditions,
                                           constraints);
  if (dim == 1)
    VectorTools::interpolate_boundary_values(dof_handler,
                                             1,
                                             boundary_conditions,
                                             constraints);
  constraints.close();

  TrilinosWrappers::SparsityPattern sp(locally_owned_dofs, mpi_communicator);
  DoFTools::make_sparsity_pattern(
    dof_handler, sp, constraints, false, this_mpi_process);
  sp.compress();

  system_matrix.reinit(sp);

  estimated_error_per_cell.reinit(triangulation.n_active_cells());

  // print out some info:
  pcout << "Number of active cells:       " << triangulation.n_active_cells()
        << std::endl;
  pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
        << std::endl;
}



template <int dim>
void
Laplace<dim>::assemble()
{
  pcout << "Assembling...";

  hp::FEValues<dim> hp_fe_values(fe,
                                 quadrature,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

  FullMatrix<double> cell_matrix;
  Vector<double>     cell_rhs;

  std::vector<types::global_dof_index> local_dof_indices;

  for (auto &cell : dof_handler.active_cell_iterators())
    if (cell->subdomain_id() == this_mpi_process)
      {
        hp_fe_values.reinit(cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
        const unsigned int  &dofs_per_cell = fe_values.dofs_per_cell;

        local_dof_indices.resize(dofs_per_cell);
        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_rhs.reinit(dofs_per_cell);

        cell_matrix = 0.;
        cell_rhs    = 0.;

        const unsigned int n_q_points =
          hp_fe_values.get_present_fe_values().n_quadrature_points;
        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = i; j < dofs_per_cell; ++j)
                  {
                    cell_matrix(i, j) += (fe_values.shape_grad(i, q_index) *
                                          fe_values.shape_grad(j, q_index)) *
                                         fe_values.JxW(q_index);
                  }

                cell_rhs(i) +=
                  force_function.value(fe_values.quadrature_point(q_index)) *
                  fe_values.shape_value(i, q_index) * fe_values.JxW(q_index);
              }
          }

        // exploit symmetry
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = i; j < dofs_per_cell; ++j)
            cell_matrix(j, i) = cell_matrix(i, j);


        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  pcout << " done." << std::endl;
}



// #define DIRECT

template <int dim>
void
Laplace<dim>::solve()
{
  pcout << "Solving...";

  SolverControl solver_control(system_rhs.size(),
                               1e-8 * system_rhs.l2_norm(),
                               /*log_history*/ false,
                               /*log_result*/ false);

  constraints.set_zero(solution);
  constraints.set_zero(system_rhs);
#ifdef DIRECT
  std::string solver_name =
    "Amesos_Superludist"; //"Amesos_Mumps" ||  "Amesos_Klu"

  TrilinosWrappers::SolverDirect::AdditionalData additional_data(false,
                                                                 solver_name);

  TrilinosWrappers::SolverDirect solver(solver_control, additional_data);

  solver.solve(system_matrix, solution, system_rhs);

  TrilinosWrappers::MPI::Vector tmp(solution);
  const double l2 = system_matrix.residual(tmp, solution, system_rhs);
  solver_control.check(1, l2);
#else
  TrilinosWrappers::SolverCG cg(solver_control);

  TrilinosWrappers::PreconditionSSOR                 preconditioner;
  TrilinosWrappers::PreconditionSSOR::AdditionalData data(1.2);
  preconditioner.initialize(system_matrix, data);

  cg.solve(system_matrix, solution, system_rhs, preconditioner);
#endif

  constraints.distribute(solution);
  solution_locally_relevant = solution;

  pcout << " done." << std::endl;
}



template <int dim>
void
Laplace<dim>::refine_grid(const unsigned int cycle)
{
  pcout << "Refining mesh..." << std::endl;

  // 3.2. Mark cells for h-refinement
  mark_h_cells();

  // 3.3. Substitute h- for p-refinement
  substitute_h_for_p();

  // prepare refinement and store number of flagged cells
  triangulation.prepare_coarsening_and_refinement();

  hp_number = {0, 0};
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->refine_flag_set())
        ++hp_number.first;
      if (cell->future_fe_index_set())
        ++hp_number.second;
    }

  // 3.4. Solution Transfer
  SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> soltrans(dof_handler);



  // copy current functions
  TrilinosWrappers::MPI::Vector solution_coarse;
  solution_coarse.reinit(locally_owned_dofs,
                         locally_relevant_dofs,
                         mpi_communicator);
  solution_coarse = solution;
  soltrans.prepare_for_coarsening_and_refinement(solution_coarse);

  // 3.5. h-refinement and p-refinement
  triangulation.execute_coarsening_and_refinement();

  // FIXME: some hp-strategies might need:
  // post_execute_coarsening_and_refinement();

  // 3.6. Setup
  setup_system();

  // 3.7. Solution Transfer finish
  soltrans.interpolate(solution);
}



template <int dim>
void
Laplace<dim>::calculate_error()
{
  L2_error     = 0.0;
  H1_error     = 0.0;
  Linfty_error = 0.0;

  hp::FEValues<dim> hp_fe_values_linf(fe,
                                      quadrature_infty,
                                      update_values | update_quadrature_points);
  hp::FEValues<dim> hp_fe_values(fe,
                                 quadrature,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

  std::vector<double>         values, exact_values;
  std::vector<double>         values_linf, exact_values_linf;
  std::vector<Tensor<1, dim>> gradients, exact_gradients;

  for (auto &cell : dof_handler.active_cell_iterators())
    if (cell->subdomain_id() == this_mpi_process)
      {
        hp_fe_values.reinit(cell);
        hp_fe_values_linf.reinit(cell);
        const FEValues<dim> &fe_values  = hp_fe_values.get_present_fe_values();
        const unsigned int   n_q_points = fe_values.n_quadrature_points;

        const FEValues<dim> &fe_values_linf =
          hp_fe_values_linf.get_present_fe_values();
        const unsigned int n_q_points_linf = fe_values_linf.n_quadrature_points;

        values_linf.resize(n_q_points_linf);
        exact_values_linf.resize(n_q_points_linf);

        values.resize(n_q_points);
        exact_values.resize(n_q_points);
        gradients.resize(n_q_points);
        exact_gradients.resize(n_q_points);

        fe_values.get_function_values(solution_locally_relevant, values);
        fe_values_linf.get_function_values(solution_locally_relevant,
                                           values_linf);
        fe_values.get_function_gradients(solution_locally_relevant, gradients);

        exact_solution.value_list(fe_values.get_quadrature_points(),
                                  exact_values);

        exact_solution.value_list(fe_values_linf.get_quadrature_points(),
                                  exact_values_linf);

        exact_solution.gradient_list(fe_values.get_quadrature_points(),
                                     exact_gradients);

        double cell_L2 = 0.0, cell_Linf = 0.0, cell_H1 = 0.0;

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double diff_values = exact_values[q_point] - values[q_point];
            const Tensor<1, dim> diff_grad =
              exact_gradients[q_point] - gradients[q_point];
            cell_L2 += diff_values * diff_values * fe_values.JxW(q_point);
            cell_H1 += (diff_grad * diff_grad) * fe_values.JxW(q_point);
          }

        for (unsigned int q_point = 0; q_point < n_q_points_linf; ++q_point)
          {
            cell_Linf = std::max(cell_Linf,
                                 std::abs(exact_values_linf[q_point] -
                                          values_linf[q_point]));
          }


        // calculate l2_norm() for cell-vectors for L2 and H1
        // and linfty_norm() for Linf:
        L2_error += cell_L2;
        H1_error += cell_H1;
        Linfty_error = std::max(Linfty_error, cell_Linf);
      } // end of loop over cells

  // finish l2_norm() / linfty_norm() calculation:
  L2_error     = sqrt(Utilities::MPI::sum(L2_error, mpi_communicator));
  H1_error     = sqrt(Utilities::MPI::sum(H1_error, mpi_communicator));
  Linfty_error = Utilities::MPI::max(Linfty_error, mpi_communicator);
}



template <int dim>
void
Laplace<dim>::output_results(int cycle)
{
  // log:
  error_table.add_value("cycle", cycle);
  error_table.add_value("cells", triangulation.n_active_cells());
  error_table.add_value("h-cells", hp_number.first);
  error_table.add_value("p-cells", hp_number.second);
  error_table.add_value("dofs", dof_handler.n_dofs());
  error_table.add_value("L2", L2_error);
  error_table.add_value("H1", H1_error);
  error_table.add_value("Linfty", Linfty_error);
  error_table.add_value("estimated", total_error);

  if (this_mpi_process == 0)
    deallog << cycle << sp << triangulation.n_active_cells() << sp
            << hp_number.first << sp << hp_number.second << sp
            << dof_handler.n_dofs() << sp << L2_error << sp << H1_error << sp
            << Linfty_error << sp << total_error << sp << std::endl;
}



template <int dim>
void
Laplace<dim>::print_errors()
{
  error_table.set_precision("L2", 3);
  error_table.set_precision("H1", 3);
  error_table.set_precision("Linfty", 3);
  error_table.set_precision("estimated", 3);
  error_table.set_scientific("L2", true);
  error_table.set_scientific("H1", true);
  error_table.set_scientific("Linfty", true);
  error_table.set_scientific("estimated", true);

  pcout << std::endl << "Error analysis:" << std::endl;
  if (this_mpi_process == 0)
    {
      error_table.write_text(std::cout);

      const std::string fname = output_name + ".gp";
      std::ofstream     output(fname, std::ios::out | std::ios::trunc);

      // use Gnuplot datablocks:
      output << "$data << EOD" << std::endl;
      error_table.write_text(output);
      output << "EOD" << std::endl << std::endl;

      output
        << "set terminal postscript eps enhanced color dashed \"Helvetica\" 22"
        << std::endl
        << "set style line 1  linetype 1 linecolor rgb \"#e41a1c\"  linewidth 2.000 pointtype 4 pointsize 2.0"
        << std::endl
        << "set style line 2  linetype 1 linecolor rgb \"#377eb8\"  linewidth 2.000 pointtype 6 pointsize 2.0"
        << std::endl
        << "set xlabel \"DoF\"" << std::endl
        << "set ylabel \"L2+H1\"" << std::endl
        << "set logscale xy" << std::endl
        << "set format x \"10^{%T}\"" << std::endl
        << "set format y \"10^{%T}\"" << std::endl
        << "set output \'" << output_name << ".eps\'" << std::endl
        << "plot \"$data\" using ($5):($6+$7) axis x1y1 with lp ls 1 title  \"error\", \\"
        << std::endl
        << "     \"$data\" using ($5):($9) axis x1y1 with lp ls 2 title     \"{/Symbol h}_{/Symbol W}\""
        << std::endl;
    }
}



template <int dim>
void
Laplace<dim>::run()
{
  // 1. Define problem
  setup_geometry();
  setup_system();

  for (unsigned int cycle = 0; cycle <= n_cycles; ++cycle)
    {
      pcout << std::endl << "Cycle " << cycle << std::endl;

      // 2. Solve Problem
      assemble();
      solve();
      calculate_error();

      estimate_error();
      total_error = estimated_error_per_cell.l2_norm();

      output_results(cycle);

      // Do refinement (Yes/No) ?
      if (cycle < n_cycles)
        {
          refine_grid(cycle);
        }
    }
  print_errors();
}
