/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Toby D. Young, Polish Academy of Sciences,
 *          Wolfgang Bangerth, Texas A&M University
 */

#include "../tests.h"

#include <deal.II/base/logstream.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/lac/arpack_solver.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <complex>

#include <fstream>
#include <iostream>
#include <algorithm>

namespace Step36
{
  using namespace dealii;


  template <int dim>
  class EigenvalueProblem
  {
  public:
    EigenvalueProblem (const std::string &prm_file);
    void run ();

  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    std::pair<unsigned int, double> solve ();
    void output_results () const;

    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;

    SparsityPattern                     sparsity_pattern;
    SparseMatrix<double>                stiffness_matrix, mass_matrix;
    std::vector<Vector<double> >        eigenfunctions;
    std::vector<std::complex<double> >  eigenvalues;

    ConstraintMatrix constraints;
  };



  template <int dim>
  EigenvalueProblem<dim>::EigenvalueProblem (const std::string &prm_file)
    :
    fe (1),
    dof_handler (triangulation)
  {
  }



  template <int dim>
  void EigenvalueProblem<dim>::make_grid_and_dofs ()
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
    triangulation.refine_global (5);
    dof_handler.distribute_dofs (fe);

    DoFTools::make_zero_boundary_constraints (dof_handler, constraints);
    constraints.close ();

    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    constraints.condense (sparsity_pattern);
    sparsity_pattern.compress();
    stiffness_matrix.reinit (sparsity_pattern);
    mass_matrix.reinit (sparsity_pattern);

    eigenfunctions.resize (8);
    for (unsigned int i=0; i<eigenfunctions.size (); ++i)
      eigenfunctions[i].reinit (dof_handler.n_dofs ());

    eigenvalues.resize (eigenfunctions.size ());
  }



  template <int dim>
  void EigenvalueProblem<dim>::assemble_system ()
  {
    QGauss<dim>   quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_stiffness_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_mass_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    FunctionParser<dim> potential;
    potential.initialize (FunctionParser<dim>::default_variable_names (),
                          "0",
                          typename FunctionParser<dim>::ConstMap());

    std::vector<double> potential_values (n_q_points);


    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active (),
    endc = dof_handler.end ();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        cell_stiffness_matrix = 0;
        cell_mass_matrix      = 0;

        potential.value_list (fe_values.get_quadrature_points(),
                              potential_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                cell_stiffness_matrix (i, j)
                += (fe_values.shape_grad (i, q_point) *
                    fe_values.shape_grad (j, q_point)
                    +
                    potential_values[q_point] *
                    fe_values.shape_value (i, q_point) *
                    fe_values.shape_value (j, q_point)
                   ) * fe_values.JxW (q_point);

                cell_mass_matrix (i, j)
                += (fe_values.shape_value (i, q_point) *
                    fe_values.shape_value (j, q_point)
                   ) * fe_values.JxW (q_point);
              }

        cell->get_dof_indices (local_dof_indices);

        constraints
        .distribute_local_to_global (cell_stiffness_matrix,
                                     local_dof_indices,
                                     stiffness_matrix);
        constraints
        .distribute_local_to_global (cell_mass_matrix,
                                     local_dof_indices,
                                     mass_matrix);
      }

    stiffness_matrix.compress (VectorOperation::add);
    mass_matrix.compress (VectorOperation::add);


    double min_spurious_eigenvalue = std::numeric_limits<double>::max(),
           max_spurious_eigenvalue = -std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
      if (constraints.is_constrained(i))
        {
          const double ev = stiffness_matrix(i,i)/mass_matrix(i,i);
          min_spurious_eigenvalue = std::min (min_spurious_eigenvalue, ev);
          max_spurious_eigenvalue = std::max (max_spurious_eigenvalue, ev);
        }
  }



  template <int dim>
  std::pair<unsigned int, double> EigenvalueProblem<dim>::solve ()
  {
    SolverControl solver_control (dof_handler.n_dofs(), 1e-10);
    SparseDirectUMFPACK inverse;
    inverse.initialize (stiffness_matrix);
    const unsigned int num_arnoldi_vectors = 2*eigenvalues.size() + 2;
    ArpackSolver::AdditionalData additional_data(num_arnoldi_vectors,
                                                 ArpackSolver::largest_real_part);
    ArpackSolver eigensolver (solver_control, additional_data);
    eigensolver.solve (stiffness_matrix,
                       mass_matrix,
                       inverse,
                       eigenvalues,
                       eigenfunctions,
                       eigenvalues.size());

    // make sure that we have eigenvectors and they are mass-orthonormal:
    // a) (A*x_i-\lambda*B*x_i).L2() == 0
    // b) x_j*B*x_i=\delta_{ij}
    {
      Vector<double> Ax(eigenfunctions[0]), Bx(eigenfunctions[0]);
      for (unsigned int i=0; i < eigenfunctions.size(); ++i)
        {
          mass_matrix.vmult(Bx,eigenfunctions[i]);

          for (unsigned int j=0; j < eigenfunctions.size(); j++)
            Assert( std::abs( eigenfunctions[j] * Bx - (i==j))< 1e-8,
                    ExcMessage("Eigenvectors " +
                               Utilities::int_to_string(i) +
                               " and " +
                               Utilities::int_to_string(j) +
                               " are not orthonormal!"
                               " failing norm is " +
                               Utilities::to_string(
                                 std::abs( eigenfunctions[j] * Bx - (i==j) )
                               )
                              ));

          stiffness_matrix.vmult(Ax,eigenfunctions[i]);
          Ax.add(-1.0*std::real(eigenvalues[i]),Bx);
          Assert (Ax.l2_norm() < 1e-8,
                  ExcMessage("Returned vector " +
                             Utilities::int_to_string(i) +
                             " is not an eigenvector!"
                             " L2 norm of the residual is " +
                             Utilities::to_string(Ax.l2_norm())
                            ));
        }
    }
    for (unsigned int i=0; i<eigenfunctions.size(); ++i)
      eigenfunctions[i] /= eigenfunctions[i].linfty_norm ();

    return std::make_pair (solver_control.last_step (),
                           solver_control.last_value ());
  }



  template <int dim>
  void EigenvalueProblem<dim>::output_results () const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);

    for (unsigned int i=0; i<eigenfunctions.size(); ++i)
      data_out.add_data_vector (eigenfunctions[i],
                                std::string("eigenfunction_") +
                                Utilities::int_to_string(i));

    Vector<double> projected_potential (dof_handler.n_dofs());
    {
      FunctionParser<dim> potential;
      potential.initialize (FunctionParser<dim>::default_variable_names (),
                            "0",
                            typename FunctionParser<dim>::ConstMap());
      VectorTools::interpolate (dof_handler, potential, projected_potential);
    }
    data_out.add_data_vector (projected_potential, "interpolated_potential");

    data_out.build_patches ();

    std::ofstream output ("eigenvectors.vtk");
    data_out.write_vtk (output);
  }


  bool my_compare(std::complex<double> a, std::complex<double> b)
  {
    return a.real() < b.real();
  }

  template <int dim>
  void EigenvalueProblem<dim>::run ()
  {
    make_grid_and_dofs ();

    assemble_system ();

    const std::pair<unsigned int, double> res = solve ();

    std::sort(eigenvalues.begin(), eigenvalues.end(), my_compare);

    for (unsigned int i = 0; i < 5 && i < eigenvalues.size(); ++i)
      deallog << "      Eigenvalue " << i
              << " : " << eigenvalues[i]
              << std::endl;
  }
}

int main (int argc, char **argv)
{

  try
    {
      using namespace dealii;
      using namespace Step36;

      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.threshold_double(1.e-10);


      EigenvalueProblem<2> problem ("");
      problem.run ();
    }

  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
