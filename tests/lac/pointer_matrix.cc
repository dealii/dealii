//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// Test vmult and Tvmult of PointerMatrix and TransposeMatrix

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/transpose_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>

int main()
{
  std::ofstream logfile("pointer_matrix/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  FullMatrix<double> A(5,5);
  unsigned int k=0;
  for (unsigned int i=0;i<A.m();++i)
    for (unsigned int j=0;j<A.n();++j)
      A(i,j) = ++k;

  PointerMatrix<FullMatrix<double>, Vector<double> > P(&A, "P");
  TransposeMatrix<FullMatrix<double>, Vector<double> > T(&A, "T");
  
  Vector<double> x(5);
  Vector<double> y(5);
  Vector<double> y2(5);
  Vector<double> diff(5);
  
  for (unsigned int i=0;i<x.size();++i)
    {
      x = 0.;
      x(i) = 1.;
      A.vmult(y,x);
      P.vmult(y2,x);
      diff = y;
      diff.add(-1., y2);
      deallog << "P.vmult:  diff " << diff.l2_norm() << std::endl;
      T.Tvmult(y2,x);
      diff = y;
      diff.add(-1., y2);
      deallog << "T.Tvmult: diff " << diff.l2_norm() << std::endl;
      
      A.Tvmult(y,x);
      P.Tvmult(y2,x);
      diff = y;
      diff.add(-1., y2);
      deallog << "P.Tvmult: diff " << diff.l2_norm() << std::endl;
      T.vmult(y2,x);
      diff = y;
      diff.add(-1., y2);
      deallog << "T.vmult:  diff " << diff.l2_norm() << std::endl;
      
      
    }
}
