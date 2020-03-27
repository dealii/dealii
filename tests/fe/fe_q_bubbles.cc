/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2018 by the deal.II authors
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
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 */

// Test whether the bubble functions can be approximated exactly

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_bubbles.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


template <int dim>
class BubbleFunction : public Function<dim>
{
public:
  BubbleFunction(unsigned int degree, unsigned int direction);

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;

  virtual Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component = 0) const;

private:
  unsigned int m_degree;
  unsigned int m_direction;
};

template <int dim>
BubbleFunction<dim>::BubbleFunction(unsigned int degree, unsigned int direction)
  : Function<dim>()
  , m_degree(degree)
  , m_direction(direction)
{}

template <int dim>
double
BubbleFunction<dim>::value(const Point<dim> &p, const unsigned int) const
{
  double return_value = 1.;
  for (unsigned int i = 0; i < dim; ++i)
    return_value *= (1 - p(i) * p(i));
  return_value *= std::pow(p(m_direction), m_degree - 1);

  return return_value;
}

template <int dim>
Tensor<1, dim>
BubbleFunction<dim>::gradient(const Point<dim> &p, const unsigned int) const
{
  Tensor<1, dim> grad;

  for (unsigned int d = 0; d < dim; ++d)
    {
      grad[d] = 1.;
      // compute grad(\prod_{i=1}^d (1-x_i^2))(p)
      for (unsigned j = 0; j < dim; ++j)
        grad[d] *= (d == j ? -2 * p(j) : (1 - p(j) * p(j)));
      // and multiply with x_i^{r-1}
      grad[d] *= std::pow(p(m_direction), m_degree - 1);
    }

  if (m_degree >= 2)
    {
      // add \prod_{i=1}^d (1-x_i^2))(p)
      double value = 1.;
      for (unsigned int j = 0; j < dim; ++j)
        value *= (1 - p(j) * p(j));
      // and multiply with grad(x_i^{r-1})
      grad[m_direction] +=
        value * (m_degree - 1) * std::pow(p(m_direction), m_degree - 2);
    }

  return grad;
}


template <int dim>
class Step3
{
public:
  Step3(FiniteElement<dim> *fe, const unsigned int degree);

  void
  run();


private:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_system(unsigned int i);
  void
  solve();
  void
  output_results(unsigned int i) const;

  Triangulation<dim>  triangulation;
  FiniteElement<dim> *fe;
  DoFHandler<dim>     dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;

  const unsigned int m_degree;
};

template <int dim>
Step3<dim>::Step3(FiniteElement<dim> *fe, const unsigned int degree)
  : fe(fe)
  , dof_handler(triangulation)
  , m_degree(degree + 1)
{}


template <int dim>
void
Step3<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(0);
}

template <int dim>
void
Step3<dim>::setup_system()
{
  dof_handler.distribute_dofs(*fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, c_sparsity);
  sparsity_pattern.copy_from(c_sparsity);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim>
void
Step3<dim>::assemble_system(unsigned int i)
{
  system_matrix = 0.;
  system_rhs    = 0.;
  QGauss<dim>   quadrature_formula(m_degree + 1);
  FEValues<dim> fe_values(*fe,
                          quadrature_formula,
                          update_values | update_gradients | update_JxW_values |
                            update_quadrature_points);

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  BubbleFunction<dim> bubble_function(m_degree - 1, i);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      cell_matrix = 0;
      cell_rhs    = 0;

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                (fe_values.shape_value(i, q_point) *
                 fe_values.shape_value(j, q_point) * fe_values.JxW(q_point));
            cell_rhs(i) +=
              (fe_values.shape_value(i, q_point) *
               bubble_function.value(fe_values.quadrature_point(q_point)) *
               fe_values.JxW(q_point));
          }

      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
}


template <int dim>
void
Step3<dim>::solve()
{
  SolverControl solver_control(1000, 1e-12);
  SolverCG<>    solver(solver_control);

  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}


template <int dim>
void
Step3<dim>::output_results(unsigned int i) const
{
  // Visualize the results
#ifdef DEBUG_Q_BUBBLES
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches(m_degree + 1);



  std::ofstream output(
    (fe->get_name() + "." + Utilities::int_to_string(i, 1) + ".vtk").c_str());
  data_out.write_vtk(output);
#endif

  Vector<float> difference_per_cell(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    BubbleFunction<dim>(m_degree - 1, i),
                                    difference_per_cell,
                                    QGauss<dim>(m_degree + 2),
                                    VectorTools::L2_norm);
  const double L2_error = difference_per_cell.l2_norm();
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    BubbleFunction<dim>(m_degree - 1, i),
                                    difference_per_cell,
                                    QGauss<dim>(m_degree + 2),
                                    VectorTools::H1_seminorm);
  const double         H1_error = difference_per_cell.l2_norm();
  const QTrapez<1>     q_trapez;
  const QIterated<dim> q_iterated(q_trapez, 5);
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    BubbleFunction<dim>(m_degree - 1, i),
                                    difference_per_cell,
                                    q_iterated,
                                    VectorTools::Linfty_norm);
  const double Linfty_error = difference_per_cell.linfty_norm();

  deallog << std::endl
          << fe->get_name() << " " << i << std::endl
          << "L2_error: " << L2_error << std::endl
          << "H1_error: " << H1_error << std::endl
          << "Linfty_error: " << Linfty_error << std::endl
          << std::endl;
}


template <int dim>
void
Step3<dim>::run()
{
  make_grid();
  setup_system();
  for (unsigned int i = 0; i < dim; ++i)
    {
      assemble_system(i);
      solve();
      output_results(i);
    }
}


int
main()
{
  initlog();
  deallog.depth_file(1);
  for (unsigned int degree = 1; degree <= 3; ++degree)
    {
      //     {
      //       FiniteElement<2> *fe = new FE_Q<2>(degree);
      //       {
      //         Step3<2> laplace_problem(fe, degree);
      //         laplace_problem.run();
      //        }
      //        delete fe;
      //     }

      {
        FiniteElement<2> *fe = new FE_Q_Bubbles<2>(degree);
        {
          Step3<2> laplace_problem(fe, degree);
          laplace_problem.run();
        }
        delete fe;
      }

      //     {
      //       FiniteElement<3> *fe = new FE_Q<3>(degree);
      //       {
      //         Step3<3> laplace_problem(fe, degree);
      //         laplace_problem.run();
      //        }
      //        delete fe;
      //     }

      {
        FiniteElement<3> *fe = new FE_Q_Bubbles<3>(degree);
        {
          Step3<3> laplace_problem(fe, degree);
          laplace_problem.run();
        }
        delete fe;
      }
    }
  return 0;
}
