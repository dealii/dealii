// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include "../tests.h"

// Check that reinit correctly set all the entries of the vector to zero

void
test()
{
  const unsigned int                          size = 100;
  LinearAlgebra::CUDAWrappers::Vector<double> a(size);
  LinearAlgebra::ReadWriteVector<double>      read_write(size);
  for (unsigned int i = 0; i < size; ++i)
    read_write[i] = i;
  a.import(read_write, VectorOperation::insert);

  a.reinit(size / 2);
  AssertThrow(a.l1_norm() == 0., ExcMessage("reinit did not zero the entry"));
}


int
main(int argc, char **argv)
{
  initlog();
  deallog.depth_console(0);

  test();

  deallog << "OK" << std::endl;

  return 0;
}
