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
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  double a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};
  double b[3][3] = {{25,31,37}, {45,57,69}, {75,96,117}};

  const unsigned int dim=3;
  Tensor<2,dim> t(a);
  Tensor<2,dim> tt;
  Tensor<2,dim> result(b);
  Assert (transpose(transpose(t)) == t, ExcInternalError());
  Assert (transpose(transpose(result)) == result, ExcInternalError());

  Vector<double> unrolled(9);

  t.unroll(unrolled);
  deallog << "unrolled:";
  for (unsigned i=0; i<9; i++)
    deallog << ' ' << unrolled(i);
  deallog << std::endl;

  deallog << "t=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
        deallog << t[i][j] << ' ';
      deallog << std::endl;
    };

  contract (tt,t,t);

  deallog << "tt=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
        deallog << tt[i][j] << ' ';
      deallog << std::endl;
    };

  if (true)
    {
      deallog.push("Cross product");
      Tensor<1,3> e1;
      Tensor<1,3> e2;
      Tensor<1,3> e3;
      e1[0] = 1.;
      e2[1] = 1.;
      e3[2] = 1.;
      Tensor<1,3> result;
      cross_product(result,e1,e2);
      deallog << '\t' << result[0]
              << '\t' << result[1]
              << '\t' << result[2] << std::endl;

      cross_product(result,e2,e3);
      deallog << '\t' << result[0]
              << '\t' << result[1] << '\t'
              << result[2] << std::endl;

      cross_product(result,e3,e1);
      deallog << '\t' << result[0]
              << '\t' << result[1]
              << '\t' << result[2] << std::endl;

      deallog.pop();
    }

  if (tt==result)
    {
      deallog << "Result OK." << std::endl;
    }
  else
    {
      deallog << "Result WRONG!" << std::endl;
    };

}
