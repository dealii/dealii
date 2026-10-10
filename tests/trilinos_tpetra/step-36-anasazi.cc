/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2009 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Authors: Toby D. Young, Polish Academy of Sciences,
 *          Wolfgang Bangerth, Texas A&M University
 */
 
 
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
 
#include <deal.II/base/index_set.h>
 
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
 
#include <deal.II/lac/anasazi_solver.h>
 
#include <fstream>
#include <iostream>
 
namespace Step36
{
  using namespace dealii;
 
 
  template <int dim>
  class EigenvalueProblem
  {
  public:
    EigenvalueProblem(const std::string &prm_file);
    void run();
 
  private:
    void         make_grid_and_dofs();
    void         assemble_system();
    unsigned int solve();
    void         output_results() const;
 
    Triangulation<dim> triangulation;
    const FE_Q<dim>    fe;
    DoFHandler<dim>    dof_handler;
 
    LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default>   stiffness_matrix, mass_matrix;
    std::vector<LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>> eigenfunctions;
    std::vector<double>                     eigenvalues;
 
    ParameterHandler parameters;
 
    AffineConstraints<double> constraints;
  };
 
 
 
  template <int dim>
  EigenvalueProblem<dim>::EigenvalueProblem(const std::string &prm_file)
    : fe(1)
    , dof_handler(triangulation)
  {
    parameters.declare_entry(
      "Global mesh refinement steps",
      "5",
      Patterns::Integer(0, 20),
      "The number of times the 1-cell coarse mesh should "
      "be refined globally for our computations.");
    parameters.declare_entry("Number of eigenvalues/eigenfunctions",
                             "5",
                             Patterns::Integer(0, 100),
                             "The number of eigenvalues/eigenfunctions "
                             "to be computed.");
    parameters.declare_entry("Potential",
                             "0",
                             Patterns::Anything(),
                             "A functional description of the potential.");
 
    parameters.parse_input(prm_file);
  }
 
 
 
  template <int dim>
  void EigenvalueProblem<dim>::make_grid_and_dofs()
  {
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(
      parameters.get_integer("Global mesh refinement steps"));
    dof_handler.distribute_dofs(fe);
 
    DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
    constraints.close();
 
    stiffness_matrix.reinit(dof_handler.n_dofs(),
                            dof_handler.n_dofs(),
                            dof_handler.max_couplings_between_dofs());
    mass_matrix.reinit(dof_handler.n_dofs(),
                       dof_handler.n_dofs(),
                       dof_handler.max_couplings_between_dofs());
 
    IndexSet eigenfunction_index_set = dof_handler.locally_owned_dofs();
    eigenfunctions.resize(
      parameters.get_integer("Number of eigenvalues/eigenfunctions"));
    for (auto &eigenfunction : eigenfunctions)
      eigenfunction.reinit(eigenfunction_index_set, MPI_COMM_WORLD);
 
    eigenvalues.resize(eigenfunctions.size());
  }
 
 
 
  template <int dim>
  void EigenvalueProblem<dim>::assemble_system()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);
 
    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
 
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();
 
    FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 
    FunctionParser<dim> potential;
    potential.initialize(FunctionParser<dim>::default_variable_names(),
                         parameters.get("Potential"),
                         typename FunctionParser<dim>::ConstMap());
 
    std::vector<double> potential_values(n_q_points);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell_stiffness_matrix = 0;
        cell_mass_matrix      = 0;
 
        potential.value_list(fe_values.get_quadrature_points(),
                             potential_values);
 
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                cell_stiffness_matrix(i, j) +=           
                  (fe_values.shape_grad(i, q_point) *    
                     fe_values.shape_grad(j, q_point)    
                   +                                     
                   potential_values[q_point] *           
                     fe_values.shape_value(i, q_point) * 
                     fe_values.shape_value(j, q_point)   
                   ) *                                   
                  fe_values.JxW(q_point);                
 
                cell_mass_matrix(i, j) +=              
                  (fe_values.shape_value(i, q_point) * 
                   fe_values.shape_value(j, q_point)   
                   ) *                                 
                  fe_values.JxW(q_point);              
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
 
 
    double min_spurious_eigenvalue = std::numeric_limits<double>::max(),
           max_spurious_eigenvalue = std::numeric_limits<double>::lowest();
 
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
      if (constraints.is_constrained(i))
        {
          const double ev         = stiffness_matrix(i, i) / mass_matrix(i, i);
          min_spurious_eigenvalue = std::min(min_spurious_eigenvalue, ev);
          max_spurious_eigenvalue = std::max(max_spurious_eigenvalue, ev);
        }
 
    std::cout << "   Spurious eigenvalues are all in the interval " << '['
              << min_spurious_eigenvalue << ',' << max_spurious_eigenvalue
              << ']' << std::endl;
  }
 
 
 
  template <int dim>
  unsigned int EigenvalueProblem<dim>::solve()
  {
    SolverControl                    solver_control(dof_handler.n_dofs(), 1e-9);
    Teuchos::ParameterList MyPL;
    MyPL.set( "Which", "SR");
    LinearAlgebra::TpetraWrappers::AnasaziSolverBlockKrylovSchur eigensolver(solver_control, MyPL);
    

    eigensolver.is_hermitian(true);

    
    eigensolver.solve(stiffness_matrix,
                      mass_matrix,
                      eigenvalues,
                      eigenfunctions,
                      eigenfunctions.size());
 
    for (auto &eigenfunction : eigenfunctions)
      eigenfunction /= eigenfunction.linfty_norm();
 
    return solver_control.last_step();
  }
 
 
 
  template <int dim>
  void EigenvalueProblem<dim>::output_results() const
  {
    DataOut<dim> data_out;
 
    data_out.attach_dof_handler(dof_handler);
 
    for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
      data_out.add_data_vector(eigenfunctions[i],
                               "eigenfunction_" + Utilities::int_to_string(i));
 
    Vector<double> projected_potential(dof_handler.n_dofs());
    {
      FunctionParser<dim> potential;
      potential.initialize(FunctionParser<dim>::default_variable_names(),
                           parameters.get("Potential"),
                           typename FunctionParser<dim>::ConstMap());
      VectorTools::interpolate(dof_handler, potential, projected_potential);
    }
    data_out.add_data_vector(projected_potential, "interpolated_potential");
 
    data_out.build_patches();
 
    std::ofstream output("eigenvectors.vtk");
    data_out.write_vtk(output);
  }
 
 
 
  template <int dim>
  void EigenvalueProblem<dim>::run()
  {
    make_grid_and_dofs();
 
    std::cout << "   Number of active cells:       "
              << triangulation.n_active_cells() << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;
 
    assemble_system();
 
    const unsigned int n_iterations = solve();
    std::cout << "   Solver converged in " << n_iterations << " iterations."
              << std::endl;
 
    output_results();
 
    std::cout << std::endl;
    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
      std::cout << "      Eigenvalue " << i << " : " << eigenvalues[i]
                << std::endl;
  }
} // namespace Step36
 
int main(int argc, char **argv)
{
  try
    {
      using namespace dealii;
      using namespace Step36;
 
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 
 
      AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
                  ExcMessage(
                    "This program can only be run in serial, use ./step-36"));
 
      EigenvalueProblem<2> problem("step-36.prm");
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
 
  std::cout << std::endl << "   Job done." << std::endl;
 
  return 0;
}