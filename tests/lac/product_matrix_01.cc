// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2014 by the deal.II authors
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


// check ProductMatrix::vmult


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>
#include <cmath>

template<typename number>
void
check()
{
  const double Adata[] =
    {
	  2, 3,
	  4, 5
    };

  const double Bdata[] =
    {
	  0, 1,
	  2, 3
    };


  FullMatrix<float> A(2,2);
  FullMatrix<double> B(2,2);

  A.fill(Adata);
  B.fill(Bdata);

  GrowingVectorMemory<Vector<number> > mem;

  ProductMatrix<Vector<number> > AB(A,B,mem);

  Vector<number> u(2);
  Vector<number> v(2);

  u(0) = 1.;
  u(1) = 2.;

  AB.vmult(v,u);

  deallog << v(0) << '\t' << v(1) << std::endl;
}


int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<double>();
  check<float>();
}
