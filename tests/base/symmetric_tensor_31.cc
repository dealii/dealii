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

// Test SymmetricTensor::memory_consumption

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/memory_consumption.h>


template <int rank, int dim>
void check ()
{
  deallog << "dim=" << dim << ", T=double "
          << MemoryConsumption::memory_consumption (SymmetricTensor<rank,dim>())
          << std::endl;
  deallog << "dim=" << dim << ", T=complex<double> "
          << MemoryConsumption::memory_consumption (SymmetricTensor<rank,dim,std::complex<double> >())
          << std::endl;
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
