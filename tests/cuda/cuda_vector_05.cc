// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

  a.import_elements(read_write_1, VectorOperation::insert);
  b.import_elements(read_write_2, VectorOperation::insert);
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
