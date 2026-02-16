// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2003 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check the output of Testing::rand(). The period of the implemented algorithm
// is approximately 31*(2^31 - 1) = 66,571,993,057. Take a reasonable subset of
// this for testing. We only print every 100,000th value to keep
// the output size reasonably small.

#include "../tests.h"

int
main()
{
  initlog();
  unsigned int random_value = 0;
  for (unsigned long long int i = 0; i < 300000000ULL; ++i)
    {
      random_value = Testing::rand();
      if (i % 100000 == 0)
        deallog << random_value << '\n';
    }
  deallog << std::flush;

  return 0;
}
