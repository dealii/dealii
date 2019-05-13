// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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
