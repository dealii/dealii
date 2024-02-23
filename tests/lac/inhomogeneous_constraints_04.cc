// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// We take two cells with two common dofs and create the system_matrix
// and the right-hand-side-vector for this system.  But we have the following
// two inhomogeneous constraints:
//       x_1 = -5, x_3 = 2.0  and  x_4 = 0.0
// And we want to test if the distribute_local_to_global function supplies the
// same system_matrix, right-hand-side-vector and solution as we obtain by using
// the MatrixTools::apply_boundary_values function.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <iostream>

#include "../tests.h"

void
test(bool use_constraint_matrix)
{
  DynamicSparsityPattern csp(8, 8);
  for (unsigned int i = 0; i < 8; ++i)
    csp.add(i, i);

  SparsityPattern sp;
  sp.copy_from(csp);
  SparseMatrix<double> mat(sp);
  Vector<double>       rhs(8);
  Vector<double>       solution(8);

  // "assemble":

  std::vector<types::global_dof_index> local_dofs1;
  for (unsigned int i = 0; i < 5; ++i)
    local_dofs1.push_back(i);

  std::vector<types::global_dof_index> local_dofs2;
  local_dofs2.push_back(1);
  for (unsigned int i = 1; i < 5; ++i)
    local_dofs2.push_back(3 + i);

  FullMatrix<double> local_mat(5, 5);
  Vector<double>     local_vec(5);
  for (unsigned int i = 0; i < 5; ++i)
    local_mat(i, i) = 2.0;

  local_vec = 1;

  if (use_constraint_matrix == true)
    {
      AffineConstraints<double> cm;

      cm.constrain_dof_to_zero(1);
      cm.set_inhomogeneity(1, -5.0);
      cm.constrain_dof_to_zero(3);
      cm.set_inhomogeneity(3, 2.0);
      cm.constrain_dof_to_zero(4);
      cm.set_inhomogeneity(4, 0.0);

      cm.close();
      cm.print(deallog.get_file_stream());

      cm.distribute_local_to_global(
        local_mat, local_vec, local_dofs1, mat, rhs, true);
      cm.distribute_local_to_global(
        local_mat, local_vec, local_dofs2, mat, rhs, true);
    }
  else
    {
      for (unsigned int i = 0; i < 5; ++i)
        {
          mat.add(local_dofs1[i], local_dofs1[i], local_mat(i, i));
          rhs(local_dofs1[i]) += local_vec(i);
        }

      for (unsigned int i = 0; i < 5; ++i)
        {
          mat.add(local_dofs2[i], local_dofs2[i], local_mat(i, i));
          rhs(local_dofs2[i]) += local_vec(i);
        }

      std::map<types::global_dof_index, double> boundary_values;
      boundary_values.insert(
        std::pair<types::global_dof_index, double>(1, -5.0));
      boundary_values.insert(
        std::pair<types::global_dof_index, double>(3, 2.0));
      boundary_values.insert(
        std::pair<types::global_dof_index, double>(4, 0.0));
      MatrixTools::apply_boundary_values(boundary_values, mat, solution, rhs);
    }

  mat.print(deallog.get_file_stream());
  rhs.print(deallog.get_file_stream());

  for (unsigned int i = 0; i < 8; ++i)
    {
      solution(i) = rhs(i) / mat(i, i);
    }

  solution.print(deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog.get_file_stream() << std::setprecision(2);

  // Use the constraints for the right-hand-side
  {
    test(true);
  }

  // Don not use the constraints for the right-hand-side
  {
    test(false);
  }
}
