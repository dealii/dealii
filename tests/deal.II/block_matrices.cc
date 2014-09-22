// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */
/* Program is based on /examples/step-3
   Purpose: compare the results when using a normal matrix and
   a block matrix
*/

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>


template <class Vector, class Matrix, class Sparsity>
class LaplaceProblem
{
public:
  LaplaceProblem (const unsigned int n_blocks);

  void run ();
  void reinit_sparsity ();
  void reinit_vectors ();

  Vector       solution;

private:
  void make_grid_and_dofs ();
  void assemble_system ();
  void solve ();

  const unsigned int n_blocks;

  Triangulation<2>     triangulation;
  FE_Q<2>              fe;
  DoFHandler<2>        dof_handler;

  Sparsity      sparsity_pattern;
  Matrix        system_matrix;

  Vector       system_rhs;
};


template <class Vector, class Matrix, class Sparsity>
LaplaceProblem<Vector,Matrix,Sparsity>::LaplaceProblem (const unsigned int n_blocks) :
  n_blocks (n_blocks),
  fe(1),
  dof_handler (triangulation)
{
  sparsity_pattern.reinit (n_blocks, n_blocks);
}



template <>
LaplaceProblem<Vector<double>,SparseMatrix<double>,SparsityPattern>::LaplaceProblem (const unsigned int n_blocks) :
  n_blocks (n_blocks),
  fe(1),
  dof_handler (triangulation)
{}



template <>
LaplaceProblem<Vector<float>,SparseMatrix<float>,SparsityPattern>::LaplaceProblem (const unsigned int n_blocks) :
  n_blocks (n_blocks),
  fe(1),
  dof_handler (triangulation)
{}



template <class Vector, class Matrix, class Sparsity>
void LaplaceProblem<Vector,Matrix,Sparsity>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (3);
  deallog << "Number of active cells: "
          << triangulation.n_active_cells()
          << std::endl;
  deallog << "Total number of cells: "
          << triangulation.n_cells()
          << std::endl;

  dof_handler.distribute_dofs (fe);

  deallog << "Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;

  reinit_sparsity ();
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);
  reinit_vectors ();
}


template <>
void LaplaceProblem<Vector<double>,SparseMatrix<double>,SparsityPattern>::reinit_sparsity ()
{
  sparsity_pattern.reinit (dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
}



template <>
void LaplaceProblem<Vector<double>,SparseMatrix<double>,SparsityPattern>::reinit_vectors ()
{
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}



template <>
void LaplaceProblem<Vector<float>,SparseMatrix<float>,SparsityPattern>::reinit_sparsity ()
{
  sparsity_pattern.reinit (dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
}



template <>
void LaplaceProblem<Vector<float>,SparseMatrix<float>,SparsityPattern>::reinit_vectors ()
{
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}



template <>
void LaplaceProblem<BlockVector<double>,BlockSparseMatrix<double>,BlockSparsityPattern>::reinit_sparsity ()
{
  switch (n_blocks)
    {
    case 2:
    {
      const types::global_dof_index n_dofs = dof_handler.n_dofs();
      const types::global_dof_index block_size[2] = { n_dofs/3, n_dofs - n_dofs/3 };

      for (unsigned int i=0; i<2; ++i)
        for (unsigned int j=0; j<2; ++j)
          sparsity_pattern.block(i,j).reinit (block_size[i], block_size[j],
                                              dof_handler.max_couplings_between_dofs());
      sparsity_pattern.collect_sizes ();

      break;
    };

    case 3:
    {
      const types::global_dof_index n_dofs = dof_handler.n_dofs();
      const types::global_dof_index block_size[3] = { n_dofs/5, n_dofs/7, n_dofs - n_dofs/5 - n_dofs/7 };

      for (unsigned int i=0; i<3; ++i)
        for (unsigned int j=0; j<3; ++j)
          sparsity_pattern.block(i,j).reinit (block_size[i], block_size[j],
                                              dof_handler.max_couplings_between_dofs());
      sparsity_pattern.collect_sizes ();

      break;
    };

    default:
      AssertThrow (false, ExcNotImplemented());
    };
}



template <>
void LaplaceProblem<BlockVector<double>,BlockSparseMatrix<double>,BlockSparsityPattern>::reinit_vectors ()
{
  switch (n_blocks)
    {
    case 2:
    {
      const types::global_dof_index n_dofs = dof_handler.n_dofs();
      const types::global_dof_index block_size_[2] = { n_dofs/3, n_dofs - n_dofs/3 };
      const std::vector<types::global_dof_index> block_size (&block_size_[0],
                                                             &block_size_[2]);

      solution.reinit (block_size);
      system_rhs.reinit (block_size);

      break;
    };

    case 3:
    {
      const types::global_dof_index n_dofs = dof_handler.n_dofs();
      const types::global_dof_index block_size_[3] = { n_dofs/5, n_dofs/7, n_dofs - n_dofs/5 - n_dofs/7 };
      const std::vector<types::global_dof_index> block_size (&block_size_[0],
                                                             &block_size_[3]);

      solution.reinit (block_size);
      system_rhs.reinit (block_size);

      break;
    };

    default:
      AssertThrow (false, ExcNotImplemented());
    };
}



template <class Vector, class Matrix, class Sparsity>
void LaplaceProblem<Vector,Matrix,Sparsity>::assemble_system ()
{
  QGauss<2>  quadrature_formula(2);
  FEValues<2> fe_values (fe, quadrature_formula,
                         UpdateFlags(update_values    |
                                     update_gradients |
                                     update_JxW_values));

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  ::Vector<double>     cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
                                      endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      cell_matrix = 0;
      cell_rhs = 0;

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

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             cell_matrix(i,j));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    };


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


template <class Vector, class Matrix, class Sparsity>
void LaplaceProblem<Vector,Matrix,Sparsity>::solve ()
{
  SolverControl           solver_control (1000, 1e-12, false, false);
  PrimitiveVectorMemory<Vector> vector_memory;
  SolverCG<Vector>        cg (solver_control, vector_memory);

  PreconditionJacobi<Matrix> preconditioner;
  preconditioner.initialize (system_matrix, 0.8);

  cg.solve (system_matrix, solution, system_rhs,
            preconditioner);
}


template <class Vector, class Matrix, class Sparsity>
void LaplaceProblem<Vector,Matrix,Sparsity>::run ()
{
  make_grid_and_dofs ();
  assemble_system ();
  solve ();

  for (unsigned int i=0; i<solution.size(); ++i)
    deallog
    //<< typeid(Vector).name ()
    //<< ' '
    //<< typeid(Matrix).name ()
    //<< '-'
        << i << ' ' << solution(i) << std::endl;
}



int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);


  // vector of solution vectors
  std::vector<std::vector<double> > solutions;

  if (true)
    {
      LaplaceProblem<Vector<double>,SparseMatrix<double>,SparsityPattern>
      laplace_problem (2);
      laplace_problem.run ();

      solutions.push_back (std::vector<double>());
      solutions.back().resize (laplace_problem.solution.size());
      for (unsigned int i=0; i<laplace_problem.solution.size(); ++i)
        solutions.back()[i] = laplace_problem.solution(i);
    };

  if (true)
    {
      LaplaceProblem<Vector<float>,SparseMatrix<float>,SparsityPattern>
      laplace_problem (3);
      laplace_problem.run ();

      solutions.push_back (std::vector<double>());
      solutions.back().resize (laplace_problem.solution.size());
      for (unsigned int i=0; i<laplace_problem.solution.size(); ++i)
        solutions.back()[i] = laplace_problem.solution(i);
    };

  if (true)
    {
      LaplaceProblem<BlockVector<double>,BlockSparseMatrix<double>,BlockSparsityPattern>
      laplace_problem (2);
      laplace_problem.run ();

      solutions.push_back (std::vector<double>());
      solutions.back().resize (laplace_problem.solution.size());
      for (unsigned int i=0; i<laplace_problem.solution.size(); ++i)
        solutions.back()[i] = laplace_problem.solution(i);
    };

  if (true)
    {
      LaplaceProblem<BlockVector<double>,BlockSparseMatrix<double>,BlockSparsityPattern>
      laplace_problem (3);
      laplace_problem.run ();

      solutions.push_back (std::vector<double>());
      solutions.back().resize (laplace_problem.solution.size());
      for (unsigned int i=0; i<laplace_problem.solution.size(); ++i)
        solutions.back()[i] = laplace_problem.solution(i);
    };

  const unsigned int n_datasets = solutions.size();
  deallog << "Checking " << n_datasets << " data sets." << std::endl;

  for (unsigned int i=1; i<n_datasets; ++i)
    Assert (solutions[i].size() == solutions[i].size(),
            ExcInternalError());

  deallog << std::setprecision(16);
  for (unsigned int i=1; i<n_datasets; ++i)
    {
      // relative accuracy. data set
      // 1 is computed using floats
      // instead of doubles, so lower
      // our requirements
      const double accuracy = (i==1 ? 1e-6 : 1e-12);

      for (unsigned int j=0; j<solutions[0].size(); ++j)
        if ( std::fabs(solutions[i][j] - solutions[0][j]) >
             accuracy*std::fabs(solutions[i][j] + solutions[0][j]))
          {
            deallog << "Discrepancy: i=" << i << ", j=" << j
                    << ", sol[i][j]=" << solutions[i][j]
                    << ", sol[0][j]=" << solutions[0][j]
                    << std::endl;
            deallog << std::flush;
            Assert (false, ExcInternalError());
          };
    };


  return 0;
}
