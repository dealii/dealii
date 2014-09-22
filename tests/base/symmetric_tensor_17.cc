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


// compute the tensor C*unit_symmetric_tensor -- I needed this tensor for an
// application, so thought I'd make a test out of it

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  SymmetricTensor<4,dim> t;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          t[i][j][k][l] = 10000*(i==j && k==l ? 1 : 0) +
                          100 * ((i==k && j==l ? 1 : 0) +
                                 (i==l && j==k ? 1 : 0));

  deallog << "t=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          deallog << i << ' ' << j << ' ' << k << ' ' << l << ": "
                  << (int)t[i][j][k][l] << std::endl;

  // multiply t by the unit symmetric tensor
  const SymmetricTensor<2,dim> t_times_1 = t * unit_symmetric_tensor<dim>();

  deallog << "t*1=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      deallog << i << ' ' << j << ": "
              << (int)(t_times_1[i][j])
              << std::endl;

  // t_times_1 should be a multiple of the
  // unit tensor, given the structure we have
  // given to it
  Assert ((t_times_1 - (dim*10000 + 2*100)*unit_symmetric_tensor<dim>()).norm()
          < 1e-14*t_times_1.norm(),
          ExcInternalError());
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
