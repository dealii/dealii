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


// compute double contraction between two rank-4 tensors

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  SymmetricTensor<4,dim> T, A, B;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          {
            // write some entries
            // into the tensors. may
            // be overwritten by
            // subsequent writes, but
            // who cares?
            A[i][j][k][l] = (i+1)*(j+1)*(l+1)*(k+1);
            B[i][j][k][l] = (i+2)*(j+3)*(l+4)*(k+5);
          }

  T = A*B;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          {
            deallog << (int)T[i][j][k][l] << std::endl;

            // calculate result by
            // hand
            double tmp = 0;
            for (unsigned int a=0; a<dim; ++a)
              for (unsigned int b=0; b<dim; ++b)
                tmp += A[i][j][a][b] * B[a][b][k][l];

            Assert (std::fabs(T[i][j][k][l] - tmp) < 1e-14*tmp,
                    ExcInternalError());
          }
}




int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
