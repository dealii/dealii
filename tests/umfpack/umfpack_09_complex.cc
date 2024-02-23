// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// This test is like the _09 test in that it uses a real-valued
// matrix, but it then solves a linear system with a complex-valued
// right hand side and corresponding complex-valued solution vector.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

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

  BlockSparsityPattern sparsity_pattern;
  sparsity_pattern.reinit(3, 3);
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      sparsity_pattern.block(i, j).reinit(
        i != 2 ? dof_handler.n_dofs() / 3 :
                 dof_handler.n_dofs() - 2 * (dof_handler.n_dofs() / 3),
        j != 2 ? dof_handler.n_dofs() / 3 :
                 dof_handler.n_dofs() - 2 * (dof_handler.n_dofs() / 3),
        dof_handler.max_couplings_between_dofs());
  sparsity_pattern.collect_sizes();
  DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  BlockSparseMatrix<double> B;
  B.reinit(sparsity_pattern);

  {
    // for some reason, we can't
    // create block matrices directly
    // in matrixtools...
    SparsityPattern xsparsity_pattern;
    xsparsity_pattern.reinit(dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, xsparsity_pattern);
    xsparsity_pattern.compress();

    SparseMatrix<double> xB;
    xB.reinit(xsparsity_pattern);

    QGauss<dim> qr(2);
    MatrixTools::create_mass_matrix(dof_handler, qr, xB);

    // scale lower left part of the matrix by
    // 1/2 and upper right part by 2 to make
    // matrix nonsymmetric
    for (SparseMatrix<double>::iterator p = xB.begin(); p != xB.end(); ++p)
      if (p->column() < p->row())
        p->value() = p->value() / 2;
      else if (p->column() > p->row())
        p->value() = p->value() * 2;

    // check that we've done it right
    for (SparseMatrix<double>::iterator p = xB.begin(); p != xB.end(); ++p)
      if (p->column() != p->row())
        AssertThrow(xB(p->row(), p->column()) != xB(p->column(), p->row()),
                    ExcInternalError());

    // now copy stuff over
    for (SparseMatrix<double>::const_iterator i = xB.begin(); i != xB.end();
         ++i)
      B.set(i->row(), i->column(), i->value());
  }


  // for a number of different solution
  // vectors, make up a matching rhs vector
  // and check what the UMFPACK solver finds
  for (unsigned int i = 0; i < 3; ++i)
    {
      BlockVector<std::complex<double>> solution(3);
      solution.block(0).reinit(dof_handler.n_dofs() / 3);
      solution.block(1).reinit(dof_handler.n_dofs() / 3);
      solution.block(2).reinit(dof_handler.n_dofs() -
                               2 * (dof_handler.n_dofs() / 3));
      solution.collect_sizes();

      // Pick a solution vector. Normalize it in such a way that we
      // get the exact same norm as for the regular _09 test. Since
      // the actual choice of rhs shouldn't matter for this test, we
      // are free to normalize as we see fit.
      for (unsigned int j = 0; j < dof_handler.n_dofs(); ++j)
        {
          solution(j) =
            std::complex<double>(1 + j + j * (i + 1) * (i + 1), // real part
                                 1 + i +
                                   i * (j + 1) * (j + 1)); // imaginary part
          solution(j) *=
            1. * (j + j * (i + 1) * (i + 1)) / std::abs(solution(j));
        }

      // Then choose as rhs for the linear system the vector
      //   b = B*solution
      // so that later the solution of the linear system Bx=b
      // is again
      //   x = solution
      BlockVector<std::complex<double>> x;
      x.reinit(solution);
      BlockVector<std::complex<double>> b;
      b.reinit(solution);

      B.vmult(b, solution);
      x = b;
      SparseDirectUMFPACK().solve(B, x);

      // Check that we really got what we expected
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
