// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


// Test PreconditionChebyshev on more complex matrix and preconditioner


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_ilu.h>
#include "../testmatrix.h"
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cmath>



int main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);


  for (unsigned int size=4; size <= 16; size *= 2)
    {
      unsigned int dim = (size-1)*(size-1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;

      // Make matrix
      FDMatrix testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.five_point_structure(structure);
      structure.compress();
      SparseMatrix<double>  A(structure);
      testproblem.five_point(A);

      PreconditionChebyshev<SparseMatrix<double>, Vector<double>, SparseILU<double> > cheby;
      PreconditionChebyshev<SparseMatrix<double>, Vector<double>, SparseILU<double> >::AdditionalData cheby_data;
      cheby_data.preconditioner.reset(new SparseILU<double>());
      cheby_data.preconditioner->initialize(A);
      cheby_data.degree = 10;
      cheby_data.smoothing_range = 40;
      cheby.initialize(A, cheby_data);

      Vector<double> v(dim);
      Vector<double> tmp1(dim), tmp2(dim);
      for (unsigned int i=0; i<3; ++i)
        {
          for (unsigned int j=0; j<dim; ++j)
            v(j) = 1. * Testing::rand()/RAND_MAX;

          A.vmult (tmp1, v);
          cheby_data.preconditioner->vmult (tmp2, tmp1);
          tmp2 -= v;
          const double ilu_residual = tmp2.l2_norm();

          A.vmult (tmp1, v);
          cheby.vmult (tmp2, tmp1);
          tmp2 -= v;
          const double cheby_residual = tmp2.l2_norm();

          deallog << "Residual step i=" << i << ":  "
                  << " ilu=" << ilu_residual
                  << ", cheby=" << cheby_residual
                  << std::endl;
        }
    }

  return 0;
}
