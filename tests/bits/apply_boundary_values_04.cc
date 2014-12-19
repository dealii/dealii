// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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



// check MatrixTools::local_apply_boundary_values with elimination of
// columns

#include "../tests.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/fe/fe_q.h>
#include <fstream>
#include <iomanip>


void test ()
{
  const unsigned int N=4;
  SparsityPattern sparsity (N,N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      sparsity.add(i,j);
  sparsity.compress ();
  SparseMatrix<double> A(sparsity), B(sparsity);
  Vector<double> b1 (N);
  Vector<double> b2 (N);

  // then fill the two matrices and vectors
  // by setting up bogus matrix entries and
  // (1) writing them into the matrix and
  // applying boundary values later on, or
  // (2) applying them right away
  std::map<types::global_dof_index,double> boundary_values;
  boundary_values[N/2] = 42;

  // then fill the matrices
  std::vector<types::global_dof_index> local_dofs (N);
  FullMatrix<double> local_matrix (N,N);
  Vector<double> local_vector (N);
  {
    for (unsigned int i=0; i<N; ++i)
      local_dofs[i] = i;

    local_matrix = 0;
    for (unsigned int i=0; i<N; ++i)
      for (unsigned int j=0; j<N; ++j)
        local_matrix(i,j) = i*N+j;
    for (unsigned int i=0; i<N; ++i)
      local_vector(i) = i+1;

    // copy local to global by ourselves
    for (unsigned int i=0; i<N; ++i)
      for (unsigned int j=0; j<N; ++j)
        A.add (local_dofs[i], local_dofs[j], local_matrix(i,j));
    for (unsigned int i=0; i<N; ++i)
      b1(local_dofs[i]) += local_vector(i);

    // or let other functions do that after
    // removing boundary values
    MatrixTools::local_apply_boundary_values (boundary_values, local_dofs,
                                              local_matrix, local_vector,
                                              true);
    for (unsigned int i=0; i<N; ++i)
      for (unsigned int j=0; j<N; ++j)
        B.add (local_dofs[i], local_dofs[j], local_matrix(i,j));
    for (unsigned int i=0; i<N; ++i)
      b2(local_dofs[i]) += local_vector(i);
  }

  // for A, remove boundary values only now.
  Vector<double> x (N);
  MatrixTools::apply_boundary_values (boundary_values, A, x, b1, true);

  deallog << "A=" << std::endl;
  A.print_formatted (deallog.get_file_stream());
  deallog << "B=" << std::endl;
  B.print_formatted (deallog.get_file_stream());

  deallog << "b1=" << std::endl;
  b1.print (deallog.get_file_stream());
  deallog << "b2=" << std::endl;
  b2.print (deallog.get_file_stream());

  // now comes the check: we subtract B from
  // A, and make sure that the result is zero
  A.add (-1., B);
  deallog << "|A|=" << A.frobenius_norm() << std::endl;
  deallog << "|B|=" << B.frobenius_norm() << std::endl;
  Assert (A.frobenius_norm() < 1e-12*B.frobenius_norm(),
          ExcInternalError());

  // similar for b1 and b2
  b1 -= b2;
  deallog << "|b1|=" << b1.l2_norm() << std::endl;
  deallog << "|b2|=" << b2.l2_norm() << std::endl;
  Assert (b1.l2_norm() < 1e-12*b2.l2_norm(),
          ExcInternalError());
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
