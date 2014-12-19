// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



const unsigned int element_degree = 2;
const unsigned int dimension = 3;

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>


using namespace dealii;


template <int dim>
class HelmholtzProblem
{
public:
  HelmholtzProblem (const FiniteElement<dim> &fe);
  void run ();

private:
  void setup_system ();
  void assemble_system ();
  void solve ();

  Triangulation<dim>                      triangulation;
  const FiniteElement<dim>               &fe;
  DoFHandler<dim>                         dof_handler;

  ConstraintMatrix                        hanging_node_constraints;

  SparsityPattern                         sparsity_pattern;
  SparseMatrix<double>                    system_matrix;
  Vector<double>                          tri_sol, tri_rhs;
  TimerOutput                             timer;
};



template <int dim>
HelmholtzProblem<dim>::HelmholtzProblem (const FiniteElement<dim> &fe) :
  fe (fe),
  dof_handler (triangulation),
  timer(std::cout, TimerOutput::summary, TimerOutput::wall_times)
{}



template <int dim>
void HelmholtzProblem<dim>::setup_system ()
{
  timer.enter_subsection("setup mesh and matrix");

  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global(6);
  dof_handler.distribute_dofs (fe);
  std::cout << "Number of active cells:          "
            << triangulation.n_active_cells()
            << std::endl
            << "Number total degrees of freedom: "
            << dof_handler.n_dofs() << std::endl;

  hanging_node_constraints.clear ();
  IndexSet locally_relevant (dof_handler.locally_owned_dofs().size());
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant);
  hanging_node_constraints.reinit (locally_relevant);
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<dim>(),
                                            hanging_node_constraints);

  {
    CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
                                         dof_handler.n_dofs(),
                                         locally_relevant);
    DoFTools::make_sparsity_pattern (dof_handler, csp,
                                     hanging_node_constraints, false);
    sparsity_pattern.copy_from (csp);
  }
  system_matrix.reinit(sparsity_pattern);
  tri_sol.reinit (dof_handler.n_dofs());
  tri_rhs.reinit (tri_sol);

  timer.leave_subsection();
}


template <int dim>
void HelmholtzProblem<dim>::assemble_system ()
{
  timer.enter_subsection("write into matrix");

  QGauss<dim>   quadrature_formula(fe.degree+1);

  const unsigned int n_q_points    = quadrature_formula.size();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  FEValues<dim>  fe_values (fe, quadrature_formula,
                            update_values   | update_gradients |
                            update_JxW_values);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      if (fe_values.get_cell_similarity() != CellSimilarity::translation)
        cell_matrix = 0;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            if (fe_values.get_cell_similarity() != CellSimilarity::translation)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += ((fe_values.shape_grad(i,q_point) *
                                      fe_values.shape_grad(j,q_point)
                                      +
                                      fe_values.shape_value(i,q_point) *
                                      fe_values.shape_value(j,q_point)) *
                                     fe_values.JxW(q_point));
          }

      cell->get_dof_indices (local_dof_indices);
      hanging_node_constraints.distribute_local_to_global (cell_matrix,
                                                           local_dof_indices,
                                                           system_matrix);
    }

  timer.leave_subsection();
}



template <int dim>
void HelmholtzProblem<dim>::solve ()
{
  for (unsigned int i=0; i<tri_rhs.size(); i++)
    if (hanging_node_constraints.is_constrained(i)==false)
      {
        const double value = (double)myrand()/RAND_MAX;
        tri_rhs(i) = value;
      }

  timer.enter_subsection("40 matrix-vector products");

  for (unsigned int i=0; i<40; ++i)
    {
      system_matrix.vmult (tri_sol, tri_rhs);
    }

  timer.leave_subsection();
}


template <int dim>
void HelmholtzProblem<dim>::run ()
{
  setup_system();
  assemble_system();
  solve();
}

int main (int argc, char **argv)
{

  try
    {
      Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);
      deallog.depth_console (0);

      FE_Q<dimension> fe(element_degree);
      HelmholtzProblem<dimension> problem(fe);
      problem.run();
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



