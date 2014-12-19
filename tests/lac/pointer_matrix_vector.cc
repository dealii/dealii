// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// Test vmult and Tvmult of PointerMatrixVector

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Vector<double> u(5);
  Vector<double> v(5);
  Vector<double> w(1);
  Vector<float>  x(5);
  Vector<float>  y(5);
  Vector<float>  z(1);


  for (unsigned int i=0; i<u.size(); ++i)
    {
      u(i) = 1 << i;
      x(i) = 1 << i;
      v(i) = 6-i;
      y(i) = 6-i;
    }

  PointerMatrixVector<double> Mu(&u);
  Mu.vmult(w,v);
  deallog << "vmult  " << w(0) << std::endl << "Tvmult";
  w(0) = 2.;
  Mu.Tvmult(v,w);
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  PointerMatrixVector<float> Mx(&x);
  Mx.vmult(z,y);
  deallog << "vmult  " << z(0) << std::endl << "Tvmult";
  z(0) = 2.;
  Mx.Tvmult(y,z);
  for (unsigned int i=0; i<y.size(); ++i)
    deallog << ' ' << y(i);
  deallog << std::endl;
}
