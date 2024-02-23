// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// make sure that the SparseMIC applied with infinite fill-in
// generates the exact inverse matrix

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_mic.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"

// TODO:[WB] find test that is less sensitive to floating point accuracy

int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);


  for (unsigned int size = 4; size <= 16; size *= 2)
    {
      unsigned int dim = (size - 1) * (size - 1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;

      // Make matrix
      FDMatrix        testproblem(size, size);
      SparsityPattern structure(dim, dim, dim);
      //      testproblem.five_point_structure(structure);
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          structure.add(i, j);
      structure.compress();
      SparseMatrix<double> A(structure);
      testproblem.five_point(A);


      for (unsigned int test = 0; test < 2; ++test)
        {
          deallog << "Test " << test << std::endl;

          // generate sparse ILU.
          //
          // for test 1, test with
          // full pattern.  for test
          // 2, test with same
          // pattern as A
          SparsityPattern mic_pattern(dim, dim, dim);
          switch (test)
            {
              case 0:
                for (unsigned int i = 0; i < dim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    mic_pattern.add(i, j);
                break;

              case 1:
                for (unsigned int i = 0; i < dim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    if (structure(i, j) != SparsityPattern::invalid_entry)
                      mic_pattern.add(i, j);
                break;

              default:
                DEAL_II_NOT_IMPLEMENTED();
            };
          mic_pattern.compress();
          SparseMIC<double>::AdditionalData data;
          data.use_this_sparsity = &mic_pattern;
          SparseMIC<double> mic;
          mic.initialize(A, data);

          // now for three test vectors v
          // determine norm of
          // (I-BA)v, where B is MIC of A.
          // since matrix is symmetric,
          // likewise test for right
          // preconditioner
          Vector<double> v(dim);
          Vector<double> tmp1(dim), tmp2(dim);
          for (unsigned int i = 0; i < 3; ++i)
            {
              for (unsigned int j = 0; j < dim; ++j)
                v(j) = random_value<double>();

              A.vmult(tmp1, v);
              mic.vmult(tmp2, tmp1);
              tmp2 -= v;
              const double left_residual = tmp2.l2_norm() / v.l2_norm();

              mic.vmult(tmp1, v);
              A.vmult(tmp2, tmp1);
              tmp2 -= v;
              const double right_residual = tmp2.l2_norm() / v.l2_norm();


              deallog << "Relative residual with test vector " << i << ":  "
                      << " left=" << left_residual
                      << ", right=" << right_residual << std::endl;
            };
        };
    };
}
