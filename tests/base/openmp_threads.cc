// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
