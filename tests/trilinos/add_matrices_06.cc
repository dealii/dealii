// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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

// check that SparseMatrix::add(SparseMatrix) and
// SparseMatrix::copy_from(SparseMatrix) do not destroy memory when called
// from matrices with the same sparsity pattern and thus, structures depending
// on these matrices such as PreconditionJacobi do not get destroyed

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <iostream>

void
test(TrilinosWrappers::SparseMatrix& m)
{
  TrilinosWrappers::SparseMatrix m2(m.m(), m.n(), 0), m3(m.m(), m.n(), 0);

  // first set a few entries one-by-one
  for(unsigned int i = 0; i < m.m(); ++i)
    for(unsigned int j = 0; j < m.n(); ++j)
      if((i + 2 * j + 1) % 3 == 0 || i == j)
        {
          m.set(i, j, i * j * .5 + .5);
          m2.set(i, j, 1.);
          m3.set(i, j, -0.1);
        }

  m.compress(VectorOperation::insert);
  m2.compress(VectorOperation::insert);
  m3.compress(VectorOperation::insert);

  deallog << "Matrix nonzeros: " << m.n_nonzero_elements() << " "
          << m2.n_nonzero_elements() << " " << m3.n_nonzero_elements()
          << std::endl;

  m.copy_from(m2);
  m.add(0.12, m3);

  Vector<double> src(m.m()), dst(m.m());
  src = 1.;

  TrilinosWrappers::PreconditionJacobi prec;
  prec.initialize(m);

  m.vmult(dst, src);
  deallog << "Vector norm: " << dst.l2_norm() << " ";

  prec.vmult(dst, src);
  deallog << dst.l2_norm() << std::endl;

  m.copy_from(m2);
  m.add(0.06, m3);

  m.vmult(dst, src);
  deallog << "Vector norm: " << dst.l2_norm() << " ";

  prec.vmult(dst, src);
  deallog << dst.l2_norm() << std::endl;

  deallog << "OK" << std::endl;
}

int
main(int argc, char** argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      {
        TrilinosWrappers::SparseMatrix m(5U, 5U, 3U);

        test(m);
      }
    }
  catch(std::exception& exc)
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
  catch(...)
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
    };
}
