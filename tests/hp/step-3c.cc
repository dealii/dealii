// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// a tiny variation of step-3b (using up to q6, instead of only up to
// q5), but this one crashed at one time before even calling the solver


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

std::ofstream logfile("output");



class LaplaceProblem
{
public:
  LaplaceProblem ();

  void run ();

private:
  void make_grid_and_dofs ();
  void assemble_system ();
  void solve ();
  void output_results () const;

  Triangulation<2>     triangulation;
  hp::FECollection<2>              fe;
  hp::DoFHandler<2>        dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  // Although we do not have h-refinement,
  // hanging nodes will inevitably appear
  // due to different polynomial degrees.
  ConstraintMatrix     hanging_node_constraints;

  Vector<double>       solution;
  Vector<double>       system_rhs;
};


LaplaceProblem::LaplaceProblem () :
  dof_handler (triangulation)
{}



void LaplaceProblem::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (4);
  deallog << "Number of active cells: "
          << triangulation.n_active_cells()
          << std::endl;
  deallog << "Total number of cells: "
          << triangulation.n_cells()
          << std::endl;

  hp::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active (),
                                          endc = dof_handler.end ();

  unsigned int cell_no = 0;
  for (; cell != endc; ++cell)
    {
      cell->set_active_fe_index (Testing::rand() % fe.size());
      deallog << "Cell " << cell << " has fe_index=" << cell->active_fe_index()
              << std::endl;
      ++cell_no;
    }

  dof_handler.distribute_dofs (fe);
  deallog << "Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  // Create sparsity pattern.
  sparsity_pattern.reinit (dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  // Create constraints which stem from
  // the different polynomial degrees on
  // the different elements.
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);

  hanging_node_constraints.print (deallog.get_file_stream());

  hanging_node_constraints.close ();
  hanging_node_constraints.condense (sparsity_pattern);

  sparsity_pattern.compress();
  system_matrix.reinit (sparsity_pattern);
}



void LaplaceProblem::assemble_system ()
{
  hp::QCollection<2>  quadrature_formula;

  for (unsigned int p = 0; p < fe.size (); ++p)
    quadrature_formula.push_back (QGauss<2> (p + 2));

  hp::FEValues<2> x_fe_values (fe, quadrature_formula,
                               update_values | update_gradients | update_JxW_values);

  const unsigned int   max_dofs_per_cell = fe.max_dofs_per_cell ();

  FullMatrix<double>   cell_matrix (max_dofs_per_cell, max_dofs_per_cell);
  Vector<double>       cell_rhs (max_dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (max_dofs_per_cell);

  hp::DoFHandler<2>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      x_fe_values.reinit (cell);

      const FEValues<2> &fe_values = x_fe_values.get_present_fe_values();

      cell_matrix = 0;
      cell_rhs = 0;

      const unsigned int dofs_per_cell = cell->get_fe ().dofs_per_cell;
      const unsigned int   n_q_points    = quadrature_formula[cell->active_fe_index ()].size();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                 fe_values.shape_grad (j, q_point) *
                                 fe_values.JxW (q_point));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                          1 *
                          fe_values.JxW (q_point));

      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             cell_matrix(i,j));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  // Include hanging nodes.
  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);

  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<2>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}



void LaplaceProblem::solve ()
{
  SolverControl           solver_control (1000, 1e-6);
  SolverCG<>              cg (solver_control);

  deallog.depth_file(0);
  cg.solve (system_matrix, solution, system_rhs,
            PreconditionIdentity());
  deallog.depth_file(3);

  const unsigned int steps = solver_control.last_step();
  if (steps >= 629 && steps <= 630)
    {
      deallog << "Solver stopped within 629 - 630 iterations"
              << std::endl;
    }
  else
    {
      deallog << "Solver stopped after " << steps << " iterations"
              << std::endl;
    }

  solution.print (deallog.get_file_stream());

  hanging_node_constraints.distribute (solution);
}



void LaplaceProblem::output_results () const
{
  DataOut<2,hp::DoFHandler<2> > data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

  data_out.write_gnuplot (deallog.get_file_stream());
}



void LaplaceProblem::run ()
{
  FE_Q<2> fe_1 (1),
       fe_2 (2),
       fe_3 (3),
       fe_4 (4),
       fe_5 (5),
       fe_6 (6);

  fe.push_back (fe_1);
  fe.push_back (fe_2);
  fe.push_back (fe_3);
  fe.push_back (fe_4);
  fe.push_back (fe_5);
  fe.push_back (fe_6);

  make_grid_and_dofs ();
  assemble_system ();
  solve ();
  output_results ();
}



int main ()
{
  logfile.precision(6);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  LaplaceProblem laplace_problem;
  laplace_problem.run ();

  return 0;
}
