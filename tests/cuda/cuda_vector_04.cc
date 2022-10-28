// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Check LinearAlgebra::CUDAWrappers::Vector::print()

#include <deal.II/base/utilities.h>

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

void
test()
{
  const unsigned int                          size = 100;
  LinearAlgebra::CUDAWrappers::Vector<double> cuda_vector(size);

  LinearAlgebra::ReadWriteVector<double> read_write_1(size);
  for (unsigned int i = 0; i < size; ++i)
    {
      read_write_1[i] = i;
    }

  cuda_vector.import(read_write_1, VectorOperation::insert);

  cuda_vector.print(deallog.get_file_stream());
}

int
main(int argc, char **argv)
{
  initlog();
  deallog.depth_console(0);

  init_cuda();

  test();

  deallog << "OK" << std::endl;

  return 0;
}
