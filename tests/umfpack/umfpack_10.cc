// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test the umfpack sparse direct solver on a mass matrix that is slightly
// modified to make it nonsymmetric.
// this is essentially the same test as _03 but for a SparseMatrixEZ

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << dim << 'd' << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  tria.refine_global(1);

  // destroy the uniformity of the matrix by
  // refining one cell
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(8 - 2 * dim);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << "Number of dofs = " << dof_handler.n_dofs() << std::endl;

  SparsityPattern sparsity_pattern;
  sparsity_pattern.reinit(dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  SparseMatrix<double> B;
  B.reinit(sparsity_pattern);

  QGauss<dim> qr(2);
  MatrixTools::create_mass_matrix(dof_handler, qr, B);

  // scale lower left part of the matrix by
  // 1/2 and upper right part by 2 to make
  // matrix nonsymmetric
  for (SparseMatrix<double>::iterator p = B.begin(); p != B.end(); ++p)
    if (p->column() < p->row())
      p->value() = p->value() / 2;
    else if (p->column() > p->row())
      p->value() = p->value() * 2;

  // check that we've done it right
  for (SparseMatrix<double>::iterator p = B.begin(); p != B.end(); ++p)
    if (p->column() != p->row())
      AssertThrow(B(p->row(), p->column()) != B(p->column(), p->row()),
                  ExcInternalError());


  SparseMatrixEZ<double> B_ez(dof_handler.n_dofs(), dof_handler.n_dofs());
  B_ez.copy_from(B);
  B.clear();

  // for a number of different solution
  // vectors, make up a matching rhs vector
  // and check what the UMFPACK solver finds
  for (unsigned int i = 0; i < 3; ++i)
    {
      Vector<double> solution(dof_handler.n_dofs());
      Vector<double> x(dof_handler.n_dofs());
      Vector<double> b(dof_handler.n_dofs());

      for (unsigned int j = 0; j < dof_handler.n_dofs(); ++j)
        solution(j) = j + j * (i + 1) * (i + 1);

      B_ez.vmult(b, solution);
      x = b;
      SparseDirectUMFPACK().solve(B_ez, x);

      x -= solution;
      deallog << "relative norm distance = " << x.l2_norm() / solution.l2_norm()
              << std::endl;
      deallog << "absolute norms = " << x.l2_norm() << ' ' << solution.l2_norm()
              << std::endl;
      Assert(x.l2_norm() / solution.l2_norm() < 1e-8, ExcInternalError());
    }
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
