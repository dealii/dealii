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


// make sure a tensor similar to the elasticity tensor for the
// isotropic case works as expected by comparing with a full tensor

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  const double lambda = 5,
               mu     = 7;
  Tensor<4,dim> ts;
  Tensor<4,dim>          ta;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      {
        ta[i][j][i][j] += mu;
        ta[i][j][j][i] += mu;
        ta[i][i][j][j] += lambda;
      }
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          ts[i][j][k][l] = ta[i][j][k][l];

  Tensor<2,dim> as, bs;
  Tensor<2,dim>          aa, ba;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      as[i][j] = aa[i][j] = (1. + (i+1)*(j+1));

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      {
        double tmp_ij = 0;
        for (unsigned int k=0; k<dim; ++k)
          for (unsigned int l=0; l<dim; ++l)
            {
              deallog << i << ' ' << j << ' ' << k << ' ' << l << ": "
                      << ta[i][j][k][l] << ' ' << ts[i][j][k][l] << ' '
                      << aa[k][l] << ' ' << as[k][l]
                      << std::endl;
              tmp_ij += ts[i][j][k][l] * as[k][l];
            }
        bs[i][j] = tmp_ij;
      }
  double_contract (ba, ta, aa);

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      {
        Assert (as[i][j] == aa[i][j], ExcInternalError());
        Assert (bs[i][j] == ba[i][j], ExcInternalError());
      }

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
