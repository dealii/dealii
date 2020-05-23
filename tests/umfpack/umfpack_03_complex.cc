// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Like the _03 test, but with complex matrix and vector. That said,
// the matrix just so happens to have zero imaginary components, even
// though we store its entries as complex numbers.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


template <int dim>
void
test(const bool transpose = false)
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

  SparseMatrix<std::complex<double>> B;
  B.reinit(sparsity_pattern);

  // Create a real sparse matrix and copy it onto B
  {
    SparseMatrix<double> BB;
    BB.reinit(sparsity_pattern);

    QGauss<dim> qr(2);
    MatrixTools::create_mass_matrix(dof_handler, qr, BB);

    // scale lower left part of the matrix by
    // 1/2 and upper right part by 2 to make
    // matrix nonsymmetric
    for (SparseMatrix<double>::iterator p = BB.begin(); p != BB.end(); ++p)
      if (p->column() < p->row())
        p->value() = p->value() / 2;
      else if (p->column() > p->row())
        p->value() = p->value() * 2;

    // check that we've done it right
    for (SparseMatrix<double>::iterator p = BB.begin(); p != BB.end(); ++p)
      if (p->column() != p->row())
        AssertThrow(B(p->row(), p->column()) != BB(p->column(), p->row()),
                    ExcInternalError());

    for (SparseMatrix<double>::iterator p = BB.begin(); p != BB.end(); ++p)
      B.set(p->row(), p->column(), {p->value(), 0});
  }


  // for a number of different solution
  // vectors, make up a matching rhs vector
  // and check what the UMFPACK solver finds
  for (unsigned int i = 0; i < 3; ++i)
    {
      Vector<std::complex<double>> solution(dof_handler.n_dofs());
      Vector<std::complex<double>> x(dof_handler.n_dofs());
      Vector<std::complex<double>> b(dof_handler.n_dofs());

      for (unsigned int j = 0; j < dof_handler.n_dofs(); ++j)
        solution(j) =
          1. * (j + j * (i + 1) * (i + 1)) *
          std::complex<double>(1. / std::sqrt(2.), 1. / std::sqrt(2.));

      if (transpose)
        B.Tvmult(b, solution);
      else
        B.vmult(b, solution);

      x = b;
      SparseDirectUMFPACK().solve(B, x, transpose);

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

  test<1>(/*transpose =*/true);
  test<2>(/*transpose =*/true);
  test<3>(/*transpose =*/true);
}
