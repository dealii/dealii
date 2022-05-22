// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2022 by the deal.II authors
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

// test the testsuite framework.

// Verify that the "threads=N" test annotation actually works.

#include "../tests.h"

int
main()
{
  std::cout << "Number of threads " << dealii::MultithreadInfo::n_threads()
            << std::endl;
  return 0;
}
