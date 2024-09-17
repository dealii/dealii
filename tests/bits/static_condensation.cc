// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// a un-hp-ified version of hp/step-7

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_trace.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <typeinfo>

#include "../tests.h"


template <int dim>
class SolutionBase
{
protected:
  static const unsigned int n_source_centers = 3;
  static const Point<dim>   source_centers[n_source_centers];
  static const double       width;
};


template <>
const Point<1>
  SolutionBase<1>::source_centers[SolutionBase<1>::n_source_centers] =
    {Point<1>(-1.0 / 3.0), Point<1>(0.0), Point<1>(+1.0 / 3.0)};

template <>
const Point<2>
  SolutionBase<2>::source_centers[SolutionBase<2>::n_source_centers] =
    {Point<2>(-0.5, +0.5), Point<2>(-0.5, -0.5), Point<2>(+0.5, -0.5)};

template <>
const Point<3>
  SolutionBase<3>::source_centers[SolutionBase<3>::n_source_centers] = {
    Point<3>(-0.5, +0.5, 0.5),
    Point<3>(-0.5, -0.5, -0.5),
    Point<3>(+0.5, -0.5, 0)};

template <int dim>
const double SolutionBase<dim>::width = 1. / 3.;



template <int dim>
class Solution : public Function<dim>, protected SolutionBase<dim>
{
public:
  Solution()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;

  virtual Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component = 0) const;
};


template <int dim>
double
Solution<dim>::value(const Point<dim> &p, const unsigned int) const
{
  double return_value = 0;
  for (unsigned int i = 0; i < this->n_source_centers; ++i)
    {
      const Tensor<1, dim> x_minus_xi = p - this->source_centers[i];
      return_value +=
        std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
    }

  return return_value;
}


template <int dim>
Tensor<1, dim>
Solution<dim>::gradient(const Point<dim> &p, const unsigned int) const
{
  Tensor<1, dim> return_value;

  for (unsigned int i = 0; i < this->n_source_centers; ++i)
    {
      const Tensor<1, dim> x_minus_xi = p - this->source_centers[i];

      return_value +=
        (-2 / (this->width * this->width) *
         std::exp(-x_minus_xi.norm_square() / (this->width * this->width)) *
         x_minus_xi);
    }

  return return_value;
}



template <int dim>
class RightHandSide : public Function<dim>, protected SolutionBase<dim>
{
public:
  RightHandSide()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;
};


template <int dim>
double
RightHandSide<dim>::value(const Point<dim> &p, const unsigned int) const
{
  double return_value = 0;
  for (unsigned int i = 0; i < this->n_source_centers; ++i)
    {
      const Tensor<1, dim> x_minus_xi = p - this->source_centers[i];

      return_value +=
        ((2 * dim -
          4 * x_minus_xi.norm_square() / (this->width * this->width)) /
         (this->width * this->width) *
         std::exp(-x_minus_xi.norm_square() / (this->width * this->width)));
      return_value +=
        std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
    }

  return return_value;
}



template <int dim>
class HelmholtzProblem
{
public:
  enum RefinementMode
  {
    global_refinement,
    adaptive_refinement
  };

  HelmholtzProblem(const unsigned int   fe_degree,
                   const RefinementMode refinement_mode);

  void
  run();

private:
  void
  setup_system();
  void
  assemble_system(const bool do_reconstruct);
  void
  solve();
  void
  refine_grid();
  void
  process_solution(const unsigned int cycle);

  Triangulation<dim>        triangulation;
  FE_Q<dim>                 fe;
  DoFHandler<dim>           dof_handler;
  AffineConstraints<double> constraints;
  SparsityPattern           sparsity_pattern;
  SparseMatrix<double>      system_matrix;
  Vector<double>            solution;
  Vector<double>            system_rhs;

  FE_TraceQ<dim>            fe_trace;
  DoFHandler<dim>           dof_handler_trace;
  AffineConstraints<double> constraints_trace;
  SparsityPattern           sparsity_pattern_trace;
  SparseMatrix<double>      system_matrix_trace;
  Vector<double>            solution_trace_full;
  Vector<double>            solution_trace;
  Vector<double>            system_rhs_trace;

  const RefinementMode refinement_mode;

  ConvergenceTable convergence_table;
};



template <int dim>
HelmholtzProblem<dim>::HelmholtzProblem(const unsigned int   fe_degree,
                                        const RefinementMode refinement_mode)
  : fe(fe_degree)
  , dof_handler(triangulation)
  , fe_trace(fe_degree)
  , dof_handler_trace(triangulation)
  , refinement_mode(refinement_mode)
{
  deallog << "Solving with Q" << fe_degree << " elements, "
          << (refinement_mode == global_refinement ? "global" : "adaptive")
          << " refinement" << std::endl
          << "==========================================="
          << (refinement_mode == global_refinement ? "" : "==") << std::endl
          << std::endl;
}



template <int dim>
void
HelmholtzProblem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Solution<dim>(),
                                           constraints);
  constraints.close();

  {
    DynamicSparsityPattern csp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
    sparsity_pattern.copy_from(csp);
  }

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());


  dof_handler_trace.distribute_dofs(fe_trace);

  constraints_trace.clear();
  DoFTools::make_hanging_node_constraints(dof_handler_trace, constraints_trace);
  VectorTools::interpolate_boundary_values(dof_handler_trace,
                                           0,
                                           Solution<dim>(),
                                           constraints_trace);
  constraints_trace.close();

  {
    DynamicSparsityPattern csp(dof_handler_trace.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_trace,
                                    csp,
                                    constraints_trace,
                                    false);
    sparsity_pattern_trace.copy_from(csp);
  }

  system_matrix_trace.reinit(sparsity_pattern_trace);

  solution_trace.reinit(dof_handler_trace.n_dofs());
  system_rhs_trace.reinit(dof_handler_trace.n_dofs());
  solution_trace_full.reinit(dof_handler.n_dofs());

  deallog << "Number of DoFs:       " << dof_handler.n_dofs() << " / "
          << dof_handler_trace.n_dofs() << std::endl;
  deallog << "Number of matrix nnz: " << sparsity_pattern.n_nonzero_elements()
          << " / " << sparsity_pattern_trace.n_nonzero_elements() << std::endl;
}



template <int dim>
void
HelmholtzProblem<dim>::assemble_system(const bool do_reconstruct)
{
  QGauss<dim>     quadrature_formula(fe.degree + 1);
  QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

  const unsigned int n_q_points      = quadrature_formula.size();
  const unsigned int n_face_q_points = face_quadrature_formula.size();

  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  FullMatrix<double> trace_matrix(fe_trace.dofs_per_cell,
                                  fe_trace.dofs_per_cell);
  Vector<double>     trace_rhs(fe_trace.dofs_per_cell);

  FullMatrix<double> eliminate_matrix(cell_matrix.m() - trace_matrix.m(),
                                      cell_matrix.m() - trace_matrix.m());
  FullMatrix<double> temp_matrix(trace_matrix.m(), eliminate_matrix.m());

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_trace(
    fe_trace.dofs_per_cell);
  Vector<double> local_trace(fe_trace.dofs_per_cell);
  Vector<double> local_interior(dofs_per_cell - fe_trace.dofs_per_cell);

  FEValues<dim> x_fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  FEFaceValues<dim> x_fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);

  const RightHandSide<dim> right_hand_side;
  std::vector<double>      rhs_values(n_q_points);

  const Solution<dim> exact_solution;

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end(), tracec = dof_handler_trace.begin_active();
  for (; cell != endc; ++cell, ++tracec)
    {
      cell_rhs = 0;

      x_fe_values.reinit(cell);
      const FEValues<dim> &fe_values = x_fe_values.get_present_fe_values();

      right_hand_side.value_list(fe_values.get_quadrature_points(), rhs_values);

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              double                sum          = 0;
              const Tensor<1, dim> *shape_grad_i = &fe_values.shape_grad(i, 0);
              const Tensor<1, dim> *shape_grad_j = &fe_values.shape_grad(j, 0);
              const double *shape_value_i        = &fe_values.shape_value(i, 0);
              const double *shape_value_j        = &fe_values.shape_value(j, 0);
              for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
                sum += (shape_grad_i[q_index] * shape_grad_j[q_index] +
                        shape_value_i[q_index] * shape_value_j[q_index]) *
                       fe_values.JxW(q_index);
              cell_matrix(i, j) = sum;
            }

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            cell_rhs(i) += (fe_values.shape_value(i, q_point) *
                            rhs_values[q_point] * fe_values.JxW(q_point));
        }

      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        if (cell->face(face)->at_boundary() &&
            (cell->face(face)->boundary_id() == 1))
          {
            x_fe_face_values.reinit(cell, face);
            const FEFaceValues<dim> &fe_face_values =
              x_fe_face_values.get_present_fe_values();

            for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
              {
                const double neumann_value =
                  (exact_solution.gradient(
                     fe_face_values.quadrature_point(q_point)) *
                   fe_face_values.normal_vector(q_point));

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  cell_rhs(i) +=
                    (neumann_value * fe_face_values.shape_value(i, q_point) *
                     fe_face_values.JxW(q_point));
              }
          }

      const unsigned int shift = fe_trace.dofs_per_cell;
      const unsigned int sizes = fe.dofs_per_cell - fe_trace.dofs_per_cell;
      for (unsigned int i = 0; i < sizes; ++i)
        for (unsigned int j = 0; j < sizes; ++j)
          eliminate_matrix(i, j) = cell_matrix(i + shift, j + shift);
      if (sizes > 0)
        eliminate_matrix.gauss_jordan();

      cell->get_dof_indices(local_dof_indices);
      if (do_reconstruct == false)
        {
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);

          for (unsigned int i = 0; i < shift; ++i)
            for (unsigned int j = 0; j < sizes; ++j)
              {
                double sum = 0;
                for (unsigned int k = 0; k < sizes; ++k)
                  sum += cell_matrix(i, k + shift) * eliminate_matrix(k, j);
                temp_matrix(i, j) = sum;
              }
          for (unsigned int i = 0; i < shift; ++i)
            for (unsigned int j = 0; j < shift; ++j)
              {
                double sum = 0;
                for (unsigned int k = 0; k < sizes; ++k)
                  sum += temp_matrix(i, k) * cell_matrix(shift + k, j);
                trace_matrix(i, j) = cell_matrix(i, j) - sum;
              }
          for (unsigned int i = 0; i < shift; ++i)
            {
              double sum = 0;
              for (unsigned int k = 0; k < sizes; ++k)
                sum += temp_matrix(i, k) * cell_rhs(shift + k);
              trace_rhs(i) = cell_rhs(i) - sum;
            }

          tracec->get_dof_indices(dof_indices_trace);
          constraints_trace.distribute_local_to_global(trace_matrix,
                                                       trace_rhs,
                                                       dof_indices_trace,
                                                       system_matrix_trace,
                                                       system_rhs_trace);
        }
      else
        {
          tracec->get_interpolated_dof_values(solution_trace, local_trace);
          for (unsigned int i = 0; i < shift; ++i)
            solution_trace_full(local_dof_indices[i]) = local_trace(i);
          for (unsigned int i = 0; i < sizes; ++i)
            {
              double sum = 0;
              for (unsigned int k = 0; k < shift; ++k)
                sum += cell_matrix(shift + i, k) * local_trace(k);
              local_interior(i) = cell_rhs(shift + i) - sum;
            }
          for (unsigned int i = 0; i < sizes; ++i)
            {
              double sum = 0;
              for (unsigned int j = 0; j < sizes; ++j)
                sum += eliminate_matrix(i, j) * local_interior(j);
              solution_trace_full(local_dof_indices[shift + i]) = sum;
            }
        }
    }
}



template <int dim>
void
HelmholtzProblem<dim>::solve()
{
  {
    SolverControl solver_control(1000, 1e-12);
    SolverCG<>    cg(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
  }
  {
    SolverControl solver_control(1000, 1e-12);
    SolverCG<>    cg(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix_trace, 1.2);

    cg.solve(system_matrix_trace,
             solution_trace,
             system_rhs_trace,
             preconditioner);

    constraints_trace.distribute(solution_trace);
  }
}



template <int dim>
void
HelmholtzProblem<dim>::refine_grid()
{
  switch (refinement_mode)
    {
      case global_refinement:
        {
          triangulation.refine_global(1);
          break;
        }

      case adaptive_refinement:
        {
          Vector<float> estimated_error_per_cell(
            triangulation.n_active_cells());

          std::map<types::boundary_id, const Function<dim> *> neumann_boundary;
          KellyErrorEstimator<dim>::estimate(dof_handler,
                                             QGauss<dim - 1>(3),
                                             neumann_boundary,
                                             solution,
                                             estimated_error_per_cell);

          GridRefinement::refine_and_coarsen_fixed_number(
            triangulation, estimated_error_per_cell, 0.3, 0.03);

          triangulation.execute_coarsening_and_refinement();

          break;
        }

      default:
        {
          DEAL_II_NOT_IMPLEMENTED();
        }
    }
}



template <int dim>
void
HelmholtzProblem<dim>::process_solution(const unsigned int cycle)
{
  {
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::L2_norm);
    const double L2_error = difference_per_cell.l2_norm();

    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::H1_seminorm);
    const double H1_error = difference_per_cell.l2_norm();

    const unsigned int n_active_cells = triangulation.n_active_cells();
    const unsigned int n_dofs         = dof_handler.n_dofs();

    convergence_table.add_value("cycle", cycle);
    convergence_table.add_value("cells", n_active_cells);
    convergence_table.add_value("dofs", n_dofs);
    convergence_table.add_value("L2", L2_error);
    convergence_table.add_value("H1", H1_error);
  }
  {
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution_trace_full,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::L2_norm);
    const double L2_error = difference_per_cell.l2_norm();

    VectorTools::integrate_difference(dof_handler,
                                      solution_trace_full,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::H1_seminorm);
    const double H1_error = difference_per_cell.l2_norm();

    const unsigned int n_active_cells = triangulation.n_active_cells();
    const unsigned int n_dofs         = dof_handler.n_dofs();

    convergence_table.add_value("L2 sc", L2_error);
    convergence_table.add_value("H1 sc", H1_error);
  }
}



template <int dim>
void
HelmholtzProblem<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 6 - dim; ++cycle)
    {
      if (cycle == 0)
        {
          GridGenerator::hyper_cube(triangulation, -1, 1);
          triangulation.refine_global(1);

          typename Triangulation<dim>::cell_iterator cell =
                                                       triangulation.begin(),
                                                     endc = triangulation.end();
          for (; cell != endc; ++cell)
            for (const unsigned int face : GeometryInfo<dim>::face_indices())
              if ((cell->face(face)->center()[0] == -1) ||
                  (cell->face(face)->center()[1] == -1))
                cell->face(face)->set_boundary_id(1);
        }
      else
        refine_grid();


      setup_system();

      assemble_system(false);
      solve();
      assemble_system(true);

      process_solution(cycle);
    }

  convergence_table.set_precision("L2", 3);
  convergence_table.set_precision("H1", 3);

  convergence_table.set_scientific("L2", true);
  convergence_table.set_scientific("H1", true);

  convergence_table.set_precision("L2 sc", 3);
  convergence_table.set_precision("H1 sc", 3);

  convergence_table.set_scientific("L2 sc", true);
  convergence_table.set_scientific("H1 sc", true);

  deallog << std::endl;
  convergence_table.write_text(deallog.get_file_stream());
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);

  deallog.attach(logfile);

  try
    {
      for (unsigned int deg = 1; deg < 5; ++deg)
        {
          HelmholtzProblem<2> helmholtz_problem_2d(
            deg, HelmholtzProblem<2>::adaptive_refinement);

          helmholtz_problem_2d.run();

          deallog << std::endl;
        }
      for (unsigned int deg = 1; deg < 4; ++deg)
        {
          HelmholtzProblem<2> helmholtz_problem_2d(
            deg, HelmholtzProblem<2>::global_refinement);

          helmholtz_problem_2d.run();

          deallog << std::endl;
        }

      {
        HelmholtzProblem<3> helmholtz_problem_3d(
          2, HelmholtzProblem<3>::adaptive_refinement);

        helmholtz_problem_3d.run();

        deallog << std::endl;
      }
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
    }

  return 0;
}


template const double SolutionBase<2>::width;
