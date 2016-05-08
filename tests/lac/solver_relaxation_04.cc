// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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


// Compare overlapping block Jacobi relaxation with different
// permutations of the blocks. All output diffs should be zero.

#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_relaxation.h>
#include <deal.II/lac/relaxation_block.h>

#include "../tests.h"
#include "testmatrix.h"
#include <cmath>
#include <fstream>
#include <iomanip>


template<typename SolverType, typename MatrixType, typename VectorType, class PRECONDITION>
double
check_solve (SolverType         &solver,
             const MatrixType   &A,
             VectorType         &u,
             VectorType         &f,
             const PRECONDITION &P)
{
  double result = 0.;
  u = 0.;
  f = 1.;
  try
    {
      solver.solve(A,u,f,P);
    }
  catch (SolverControl::NoConvergence &e)
    {
      result = e.last_residual;
    }
  return result;
}

int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  //  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  SolverControl control(10, 1.e-3);
  SolverRelaxation<> relax(control);

  // Solve non-symmetric laplace with five-point FD
  for (unsigned int size=33; size <= 33; size *= 3)
    {
      unsigned int dim = (size-1)*(size-1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;

      // Make matrix
      FDMatrix testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.five_point_structure(structure);
      structure.compress();
      SparseMatrix<double>  A(structure);
      testproblem.five_point(A,true);

      for (unsigned int blocksize = 4; blocksize < 32; blocksize <<= 1)
        {
          deallog << "Block size " << blocksize << std::endl;

          const unsigned int n_blocks = dim/blocksize;
          RelaxationBlock<SparseMatrix<double>,double>::AdditionalData relax_data(0.7);
          RelaxationBlock<SparseMatrix<double>,double>::AdditionalData relax_data_reorder(0.7);

          relax_data.block_list.reinit(n_blocks, dim, blocksize+2);
          relax_data_reorder.block_list.reinit(n_blocks, dim, blocksize+2);
          for (unsigned int block=0; block<n_blocks; ++block)
            {
              for (int i=-1; i<(int)blocksize+1; ++i)
                if ( (int)(i+block*blocksize)>-1 && (i+block*blocksize)<dim )
                  {
                    relax_data.block_list.add(block, i+block*blocksize);
                    relax_data_reorder.block_list.add(block, i+block*blocksize);
                  }
            }
          relax_data.block_list.compress();
          relax_data_reorder.block_list.compress();

          RelaxationBlockJacobi<SparseMatrix<double>,double> relax_jacobi;
          relax_jacobi.initialize(A, relax_data);

          // reverse the order of the blocks
          relax_data_reorder.order.resize(1);
          relax_data_reorder.order[0].resize(n_blocks);
          for (unsigned int i=0; i<n_blocks; ++i)
            relax_data_reorder.order[0][i] = n_blocks-1-i;

          RelaxationBlockJacobi<SparseMatrix<double>,double> relax_jacobi_reorder;
          relax_jacobi_reorder.initialize(A, relax_data_reorder);

          Vector<double>  f(dim);
          Vector<double>  u(dim);
          Vector<double> res(dim);

          f = 1.;
          u = 1.;

          try
            {
              double r1, r2;

              deallog.push("Jacobi");
              r1 = check_solve(relax,A,u,f,relax_jacobi);
              r2 = check_solve(relax,A,u,f,relax_jacobi_reorder);
              deallog << "diff " << std::fabs(r1-r2)/r1 << std::endl;
              deallog.pop();

            }
          catch (std::exception &e)
            {
              std::cerr << "Exception: " << e.what() << std::endl;
            }
        }
    }
}