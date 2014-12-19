// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <fstream>
#include <iomanip>
#include <cmath>

template<typename number>
void check_sparse_product(const SparseMatrix<number> &m1, SparseMatrix<number> &m2)
{
  Vector<double> v(m2.n());
  Vector<double> w(m1.m());
  deallog << "Sizes\t" << v.size() << '\t' << w.size() << std::endl;

  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = i+1.;

  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  GrowingVectorMemory<Vector<double> > mem;

  ProductSparseMatrix<number, number> product(m1, m2, mem);

  product.vmult(w,v);
  for (unsigned int i=0; i<w.size(); ++i)
    deallog << ' ' << w(i);
  deallog << std::endl;

  product.vmult_add(w,v);
  for (unsigned int i=0; i<w.size(); ++i)
    deallog << ' ' << w(i);
  deallog << std::endl;

  product.Tvmult(v,w);
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  product.Tvmult_add(v,w);
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;
}


int main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  SparsityPattern sparsity1(2,3,3);
  SparsityPattern sparsity2(3,4,4);

  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<3; ++j)
      {
        sparsity1.add(i,j);
        sparsity2.add(j,i);
        sparsity2.add(j,2+i);
      }
  sparsity1.compress();
  sparsity2.compress();

  SparseMatrix<double> m1(sparsity1);
  SparseMatrix<double> m2(sparsity2);

  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<3; ++j)
      {
        m1.set(i, j, 1.*i-j);
        m2.set(j, i, 1.*i-j);
        m2.set(j, 2+i, 1.*j-i);
      }
  check_sparse_product(m1, m2);
}
