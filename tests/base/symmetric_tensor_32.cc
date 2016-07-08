// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

// Test SymmetricTensor::norm() for complex-valued tensors

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>


template <int rank, int dim>
void check ()
{
  // build a regular tensor
  SymmetricTensor<rank,dim> t;

  // build one in which all numbers are the same but purely imaginary
  SymmetricTensor<rank,dim,std::complex<double> > ti;

  // build one in which all numbers have both real and imaginary components
  SymmetricTensor<rank,dim,std::complex<double> > tc;

  for (unsigned int i=0; i<t.n_independent_components; ++i)
    {
      t.access_raw_entry(i)  = 1.0 * (i+1);
      ti.access_raw_entry(i) = std::complex<double>(0,1.0*(i+1));
      tc.access_raw_entry(i) = std::complex<double>(1.0*(i+1),1.0*(i+1));
    }

  deallog << t.norm() << ' '
          << ti.norm() << ' '
          << tc.norm() << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  deallog << "check rank 2 tensors" << std::endl;
  check<2,1>();
  check<2,2>();
  check<2,3>();

  deallog << "check rank 4 tensors" << std::endl;
  check<4,1>();
  check<4,2>();
  check<4,3>();

  deallog << "OK" << std::endl;
}
