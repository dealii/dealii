/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2020 by the deal.II authors
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
 * Authors: Toby D. Young, Polish Academy of Sciences,
 *          Wolfgang Bangerth, Texas A&M University
 *          Joscha Gedicke, U Heidelberg
 *
 * This file tests the non-symmetric interface to ARPACK for an
 advection-diffussion
 * operator. The advection is chosen in such a way that we compute complex
 eigenvalues.
 * The most critical case when we cut a complex conjugate pair is tested, i.e.
 the
 * last eigenvalue is complex but the conjugate pair is not included in the
 range
 * of eigenvalues we asked for. This is a particular case in the nonsymmetric
 arpack
 * interface that needs to be taken care of.
 *
 * We test that the computed vectors are eigenvectors and mass-normal, i.e.
 * a) (A*x_i-\lambda*B*x_i).L2() == 0
 * b) x_i*B*x_i = 1
 *
 */

#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/arpack_solver.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <complex>
#include <iostream>

#include "../tests.h"



template <int dim>
class EigenvalueProblem
{
public:
  EigenvalueProblem(unsigned int n_eigenvalues);
  void
  run();

private:
  void
  make_grid_and_dofs();
  void
  assemble_system();
  void
  solve();

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dof_handler;

  SparsityPattern                   sparsity_pattern;
  SparseMatrix<double>              stiffness_matrix, mass_matrix;
  std::vector<Vector<double>>       arpack_vectors;
  std::vector<Vector<double>>       eigenvectors;
  std::vector<std::complex<double>> eigenvalues;

  AffineConstraints<double> constraints;

  unsigned int n_eigenvalues;
};



template <int dim>
EigenvalueProblem<dim>::EigenvalueProblem(unsigned int n_eigenvalues)
  : fe(1)
  , dof_handler(triangulation)
  , n_eigenvalues(n_eigenvalues)
{}



template <int dim>
void
EigenvalueProblem<dim>::make_grid_and_dofs()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(5);
  dof_handler.distribute_dofs(fe);

  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  sparsity_pattern.reinit(dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
  constraints.condense(sparsity_pattern);
  sparsity_pattern.compress();
  stiffness_matrix.reinit(sparsity_pattern);
  mass_matrix.reinit(sparsity_pattern);

  eigenvalues.resize(n_eigenvalues);

  arpack_vectors.resize(n_eigenvalues + 1);
  for (unsigned int i = 0; i < arpack_vectors.size(); ++i)
    arpack_vectors[i].reinit(dof_handler.n_dofs());

  eigenvectors.resize(2 * n_eigenvalues);
  for (unsigned int i = 0; i < eigenvectors.size(); ++i)
    eigenvectors[i].reinit(dof_handler.n_dofs());
}



template <int dim>
void
EigenvalueProblem<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(2);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      cell_stiffness_matrix = 0;
      cell_mass_matrix      = 0;

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          const Point<dim> cur_point = fe_values.quadrature_point(q_point);
          Tensor<1, dim>   advection;
          advection[0] = 10.;
          advection[1] = 10. * cur_point[0];
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                cell_stiffness_matrix(i, j) +=
                  (fe_values.shape_grad(i, q_point) *
                     fe_values.shape_grad(j, q_point) +
                   (advection * fe_values.shape_grad(i, q_point)) *
                     fe_values.shape_value(j, q_point)) *
                  fe_values.JxW(q_point);

                cell_mass_matrix(i, j) += (fe_values.shape_value(i, q_point) *
                                           fe_values.shape_value(j, q_point)) *
                                          fe_values.JxW(q_point);
              }
        }

      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(cell_stiffness_matrix,
                                             local_dof_indices,
                                             stiffness_matrix);
      constraints.distribute_local_to_global(cell_mass_matrix,
                                             local_dof_indices,
                                             mass_matrix);
    }

  stiffness_matrix.compress(VectorOperation::add);
  mass_matrix.compress(VectorOperation::add);
}



template <int dim>
void
EigenvalueProblem<dim>::solve()
{
  SolverControl       solver_control(dof_handler.n_dofs(), 1e-10);
  SparseDirectUMFPACK inverse;
  inverse.initialize(stiffness_matrix);
  const unsigned int           num_arnoldi_vectors = 2 * eigenvalues.size() + 2;
  ArpackSolver::AdditionalData additional_data(num_arnoldi_vectors,
                                               ArpackSolver::largest_real_part);
  ArpackSolver                 eigensolver(solver_control, additional_data);
  arpack_vectors[0] = 1.;
  eigensolver.set_initial_vector(arpack_vectors[0]);
  eigensolver.solve(stiffness_matrix,
                    mass_matrix,
                    inverse,
                    eigenvalues,
                    arpack_vectors,
                    eigenvalues.size());

  // extract real and complex components of eigenvectors
  for (unsigned int i = 0; i < n_eigenvalues; ++i)
    {
      eigenvectors[i] = arpack_vectors[i];
      if (eigenvalues[i].imag() != 0.)
        {
          eigenvectors[i + eigenvalues.size()] = arpack_vectors[i + 1];
          if (i + 1 < eigenvalues.size())
            {
              eigenvectors[i + 1]                      = arpack_vectors[i];
              eigenvectors[i + 1 + eigenvalues.size()] = arpack_vectors[i + 1];
              eigenvectors[i + 1 + eigenvalues.size()] *= -1;
              ++i;
            }
        }
    }

  // make sure that we have eigenvectors and they are mass-normal:
  // a) (A*x_i-\lambda*B*x_i).L2() == 0
  // b) x_i*B*x_i = 1
  {
    Vector<double> Ax(eigenvectors[0]), Bx(eigenvectors[0]);
    Vector<double> Ay(eigenvectors[0]), By(eigenvectors[0]);
    for (unsigned int i = 0; i < n_eigenvalues; ++i)
      {
        stiffness_matrix.vmult(Ax, eigenvectors[i]);
        stiffness_matrix.vmult(Ay, eigenvectors[i + n_eigenvalues]);
        mass_matrix.vmult(Bx, eigenvectors[i]);
        mass_matrix.vmult(By, eigenvectors[i + n_eigenvalues]);

        Ax.add(-1.0 * std::real(eigenvalues[i]), Bx);
        Ax.add(std::imag(eigenvalues[i]), By);
        Ay.add(-1.0 * std::real(eigenvalues[i]), By);
        Ay.add(-1.0 * std::imag(eigenvalues[i]), Bx);
        Vector<double> tmpx(Ax), tmpy(Ay);
        tmpx.scale(Ax);
        tmpy.scale(Ay);
        tmpx += tmpy;
        if (std::sqrt(tmpx.l1_norm()) > 1e-8)
          deallog << "Returned vector " << i << " is not an eigenvector!"
                  << " L2 norm of the residual is " << std::sqrt(tmpx.l1_norm())
                  << std::endl;

        const double tmp =
          std::abs(eigenvectors[i] * Bx + eigenvectors[i + n_eigenvalues] * By -
                   1.) +
          std::abs(eigenvectors[i + n_eigenvalues] * Bx - eigenvectors[i] * By);
        if (tmp > 1e-8)
          deallog << "Eigenvector " << i << " is not normal! failing norm is"
                  << tmp << std::endl;
      }
  }
}

bool
my_compare(std::complex<double> a, std::complex<double> b)
{
  if (a.imag() == 0.)
    return a.imag() < a.imag();
  else
    return a.real() < b.real();
}

template <int dim>
void
EigenvalueProblem<dim>::run()
{
  make_grid_and_dofs();

  assemble_system();

  solve();

  std::sort(eigenvalues.begin(), eigenvalues.end(), my_compare);

  for (unsigned int i = 0; i < n_eigenvalues; ++i)
    deallog << "      Eigenvalue " << i << " : " << eigenvalues[i] << std::endl;
}


int
main(int argc, char **argv)
{
  try
    {
      initlog();


      EigenvalueProblem<2> problem(4);
      problem.run();
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
