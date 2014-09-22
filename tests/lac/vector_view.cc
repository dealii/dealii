// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>
#include <cmath>
#include <fstream>
#include <iomanip>

const unsigned int N=10;
unsigned int check_point = 0;

template <typename number>
void print (const Vector<number> &v)
{
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << v(i) << '\t';
  deallog << std::endl;
}

template<typename T>
void fill( T &a)
{
  for (unsigned int i=0; i<a.size(); ++i)
    a(i) = i;
}

int main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Vector<double>  v1(N);
  fill(v1);

  deallog << "Vector" << std::endl;
  print(v1);

  VectorView<double> v2(N, v1.begin() );
  deallog << "Vector View" << std::endl;
  print(v2);

  v2(4) = 0;
  deallog << "Modified element 4" << std::endl;
  deallog << "Vector" << std::endl;
  print(v1);

  deallog << "Vector View" << std::endl;
  print(v2);

  // Const vector.
  const Vector<double> v3(v1);

  deallog << "const Vector" << std::endl;
  print(v3);

  VectorView<double> v4(N, v3.begin());
  deallog << "const Vector View" << std::endl;
  print(v4);

  v4.reinit(N, v1.begin());
  v4.reinit(N, v3.begin());
}


