//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include "../tests.h"
#include <base/logstream.h>
#include <lac/matrix_lib.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>

#include <fstream>
#include <iostream>
#include <cmath>

template<typename number>
void check_sparse_product(const SparseMatrix<number>& m1, SparseMatrix<number>& m2)
{
  Vector<double> v(m2.n());
  Vector<double> w(m1.m());
  deallog << "Sizes\t" << v.size() << '\t' << w.size() << std::endl;
  
  for (unsigned int i=0;i<v.size();++i)
    v(i) = i+1.;
  
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;
  
  GrowingVectorMemory<Vector<double> > mem;
  
  ProductSparseMatrix<number, number> product(m1, m2, mem);
  
  product.vmult(w,v);
  for (unsigned int i=0;i<w.size();++i)
    deallog << ' ' << w(i);
  deallog << std::endl;
  
  product.vmult_add(w,v);
  for (unsigned int i=0;i<w.size();++i)
    deallog << ' ' << w(i);
  deallog << std::endl;
  
  product.Tvmult(v,w);
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;
  
  product.Tvmult_add(v,w);
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v(i);
  deallog << std::endl;
}


int main()
{
  std::ofstream logfile("matrix_lib.output");
  logfile.setf(std::ios::fixed);
  logfile.precision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  SparsityPattern sparsity1(2,3,3);
  SparsityPattern sparsity2(3,4,4);

  for (unsigned int i=0;i<2;++i)
    for (unsigned int j=0;j<3;++j)
      {
	sparsity1.add(i,j);
	sparsity2.add(j,i);
	sparsity2.add(j,2+i);
      }
  sparsity1.compress();
  sparsity2.compress();

  SparseMatrix<double> m1(sparsity1);
  SparseMatrix<double> m2(sparsity2);

  for (unsigned int i=0;i<2;++i)
    for (unsigned int j=0;j<3;++j)
      {
	m1.set(i, j, 1.*i-j);
	m2.set(j, i, 1.*i-j);
	m2.set(j, 2+i, 1.*j-i);
      }
  check_sparse_product(m1, m2);
}
