// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


#include <deal.II/base/utilities.h>

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "../tests.h"

// There was a bug where add_and_dot would give the wrong result if the size of
// the vectors was greater than BLOCK_SIZE*CHUNK_SIZE.

void
test()
{
  const unsigned int                          size = 4100;
  LinearAlgebra::CUDAWrappers::Vector<double> a(size);
  LinearAlgebra::CUDAWrappers::Vector<double> b(size);
  Vector<double>                              a_host(size);
  Vector<double>                              b_host(size);

  LinearAlgebra::ReadWriteVector<double> read_write_1(size);
  LinearAlgebra::ReadWriteVector<double> read_write_2(size);

  for (unsigned int i = 0; i < size; ++i)
    {
      read_write_1[i] = i;
      read_write_2[i] = 5. + i;
      a_host[i]       = i;
      b_host[i]       = 5. + i;
    }

  a.import(read_write_1, VectorOperation::insert);
  b.import(read_write_2, VectorOperation::insert);
  AssertThrow(a.add_and_dot(2., a, b) == a_host.add_and_dot(2., a_host, b_host),
              ExcMessage("Problem in add_and_dot"));
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
