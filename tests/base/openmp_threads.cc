// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Ensure that, when executing test suite programs, the environment variable
// OMP_NUM_THREADS is set to 2.

#include "../tests.h"


int
main()
{
  initlog();

  // This might not work on Windows; always pass if we don't have getenv
  // available
#ifdef _POSIX_C_SOURCE
  const char *omp_env_variable = getenv("OMP_NUM_THREADS");
  AssertThrow(omp_env_variable != nullptr,
              ExcMessage("The environment variable OMP_NUM_THREADS should be "
                         "defined when running the test suite."));
  if (std::string(omp_env_variable) == "2")
    {
      deallog << "OK" << std::endl;
    }
#else
  deallog << "OK" << std::endl;
#endif
}
