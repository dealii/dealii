// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
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



// check PETScWrappers::MatrixBase::clear_rows () with used second argument

#include "../tests.h"
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>


void test (PETScWrappers::MatrixBase &m)
{
  Assert (m.m() != 0, ExcInternalError());
  Assert (m.n() != 0, ExcInternalError());

  typedef PETScWrappers::MatrixBase::size_type size_type;

  // build a tri-diagonal pattern
  double norm_sqr = 0;
  unsigned int nnz = 0;
  const size_type N = m.m();
  for (size_type i=0; i<N; ++i)
    {
      if (i>=5)
        {
          const double s = Testing::rand();
          m.add (i,i-5, s);
          norm_sqr += s*s;
          ++nnz;
        }

      if (i<N-5)
        {
          const double s = Testing::rand();
          m.add (i,i+5, s);
          norm_sqr += s*s;
          ++nnz;
        }

      const double s = Testing::rand();
      m.add (i,i,s);
      norm_sqr += s*s;
      ++nnz;
    }
  m.compress (VectorOperation::add);

  deallog << m.frobenius_norm() << ' ' << std::sqrt (norm_sqr)
          << std::endl;
  deallog << m.n_nonzero_elements() << ' ' << nnz << std::endl;

  Assert (std::fabs (m.frobenius_norm() - std::sqrt(norm_sqr))
          < std::fabs (std::sqrt (norm_sqr)),
          ExcInternalError());
  Assert (m.n_nonzero_elements()-nnz == 0, ExcInternalError());

  // now remove the entries of rows N/2 and
  // N/3. set diagonal entries to rnd
  const double rnd = Testing::rand();
  for (size_type i=0; i<N; ++i)
    {
      const double s = m.el(N/2,i);
      norm_sqr -= s*s;
    }
  for (size_type i=0; i<N; ++i)
    {
      const double s = m.el(N/3,i);
      norm_sqr -= s*s;
    }
  norm_sqr += 2*rnd*rnd;

  const size_type rows[2] = { N/3, N/2 };
  m.clear_rows (std::vector<size_type>(&rows[0], &rows[2]), rnd);

  deallog << m.frobenius_norm() << ' ' << std::sqrt (norm_sqr)
          << std::endl;
  deallog << m.n_nonzero_elements() << ' ' << nnz << std::endl;

  Assert (std::fabs (m.frobenius_norm() - std::sqrt(norm_sqr))
          < std::fabs (std::sqrt (norm_sqr)),
          ExcInternalError());

  // make sure that zeroing out rows does at
  // least not add new nonzero entries (it
  // may remove some, though)
  Assert (m.n_nonzero_elements() <= nnz, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        PETScWrappers::SparseMatrix v (14,14,3);
        test (v);
      }

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
    };
}
