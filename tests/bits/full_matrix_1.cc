//----------------------------  full_matrix_1.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix_1.cc  ---------------------------


// FullMatrix::copy_from could not be compiled if we copied from a
// sparse matrix. make sure this now works

#include "../tests.h"
#include <base/logstream.h>
#include <lac/sparse_matrix.h>
#include <lac/full_matrix.h>
#include <fstream>


int main () 
{
  std::ofstream logfile("full_matrix_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  const unsigned int N = 4;
  FullMatrix<double> f(N,N);

  SparsityPattern s (N,N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      s.add (i,j);
  s.compress ();

  SparseMatrix<double> sm(s);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      sm.set (i,j,i*j);

  f.copy_from (sm);

  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      {
	deallog << i << ' ' << j << ' ' << f(i,j) << std::endl;
	Assert (f(i,j) == sm(i,j), ExcInternalError());
      }

  deallog << "OK" << std::endl;
}
