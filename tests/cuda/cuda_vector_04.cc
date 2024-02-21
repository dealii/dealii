// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

  cuda_vector.import_elements(read_write_1, VectorOperation::insert);

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
