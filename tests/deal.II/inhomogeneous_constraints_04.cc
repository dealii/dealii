// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// We take two cells with two common dofs and create the system_matrix
// and the right-hand-side-vector for this system.  But we have the following
// two inhomogeneous constraints:
//       x_1 = -5, x_3 = 2.0  and  x_4 = 0.0
// And we want to test if the distribute_local_to_global function supplies the same
// system_matrix, right-hand-side-vector and solution as we obtain by using the
// MatrixTools::apply_boundary_values function.

#include "../tests.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <fstream>
#include <iostream>
#include <complex>

std::ofstream logfile("output");

using namespace dealii;

void test(bool use_constraint_matrix)
{
  CompressedSimpleSparsityPattern csp(8,8);
  for (unsigned int i=0; i<8; ++i)
    csp.add(i,i);

  SparsityPattern sp;
  sp.copy_from(csp);
  SparseMatrix<double> mat(sp);
  Vector<double> rhs(8);
  Vector<double> solution(8);

  // "assemble":

  std::vector<types::global_dof_index> local_dofs1;
  for (unsigned int i=0; i<5; ++i)
    local_dofs1.push_back(i);

  std::vector<types::global_dof_index> local_dofs2;
  local_dofs2.push_back(1);
  for (unsigned int i=1; i<5; ++i)
    local_dofs2.push_back(3+i);

  FullMatrix<double> local_mat(5,5);
  Vector<double> local_vec(5);
  for (unsigned int i=0; i<5; ++i)
    local_mat(i,i)=2.0;

  local_vec = 1;

  if (use_constraint_matrix == true)
    {
      ConstraintMatrix cm;

      cm.add_line(1);
      cm.set_inhomogeneity(1, -5.0);
      cm.add_line(3);
      cm.set_inhomogeneity(3, 2.0);
      cm.add_line(4);
      cm.set_inhomogeneity(4, 0.0);

      cm.close();
      cm.print(logfile);

      cm.distribute_local_to_global(local_mat, local_vec, local_dofs1, mat, rhs, true);
      cm.distribute_local_to_global(local_mat, local_vec, local_dofs2, mat, rhs, true);
    }
  else
    {
      for (unsigned int i=0; i<5; i++)
        {
          mat.add (local_dofs1[i], local_dofs1[i], local_mat (i,i));
          rhs (local_dofs1[i]) += local_vec (i);
        }

      for (unsigned int i=0; i<5; i++)
        {
          mat.add (local_dofs2[i], local_dofs2[i], local_mat (i,i));
          rhs (local_dofs2[i]) += local_vec (i);
        }

      std::map<types::global_dof_index,double> boundary_values;
      boundary_values.insert (std::pair<types::global_dof_index,double>(1, -5.0));
      boundary_values.insert (std::pair<types::global_dof_index,double>(3, 2.0));
      boundary_values.insert (std::pair<types::global_dof_index,double>(4, 0.0));
      MatrixTools::apply_boundary_values (boundary_values,
                                          mat,
                                          solution,
                                          rhs);
    }

  mat.print(logfile);
  rhs.print(logfile);

  for (unsigned int i=0; i<8; ++i)
    {
      solution(i)=rhs(i)/mat(i,i);
    }

  solution.print(logfile);
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // Use the constraints for the right-hand-side
  {
    test(true);
  }

  // Don not use the constraints for the right-hand-side
  {
    test(false);
  }
}
