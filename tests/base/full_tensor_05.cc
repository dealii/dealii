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


// make sure the tensor t_ijkl=delta_ik delta_jl
// actually maps a rank-2 tensor onto twice itself

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  Tensor<4,dim> t;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          t[i][j][k][l] = ((i==k) && (j==l) ? 1 : 0);

  Tensor<2,dim> a, b;
  a[0][0] = 1;
  a[1][1] = 2;
  a[0][1] = 3;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      {
        double tmp_ij = 0;
        for (unsigned int k=0; k<dim; ++k)
          for (unsigned int l=0; l<dim; ++l)
            {
              deallog << i << ' ' << j << ' ' << k << ' ' << l << ": "
                      << t[i][j][k][l] << ' ' << a[k][l]
                      << std::endl;
              tmp_ij += t[i][j][k][l] * a[k][l];
            }
        b[i][j] = tmp_ij;
      }

  Assert (a == b, ExcInternalError());

  // try the same thing with scaled
  // tensors etc
  t *= 2;
  b.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      {
        double tmp_ij = 0;
        for (unsigned int k=0; k<dim; ++k)
          for (unsigned int l=0; l<dim; ++l)
            tmp_ij += t[i][j][k][l] * a[k][l];
        b[i][j] = tmp_ij;
      }

  Assert (a == b/2, ExcInternalError());
}




int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
