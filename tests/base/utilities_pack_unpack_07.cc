// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


// tests the influence of different compression options for
// Utilities::pack/unpack
// (based upon "utilities_pack_unpack_06")


#include <deal.II/base/point.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <random>

#include "../tests.h"



template <int N>
void
check(const double (&array)[N])
{
  // ----- PACK -----
  std::vector<char> array_compressed, array_uncompressed;

  TimerOutput computing_timer(std::cout,
                              TimerOutput::never,
                              TimerOutput::wall_times);

#ifdef DEAL_II_WITH_ZLIB
  // default option should work for compression
  {
    TimerOutput::Scope timer_section(computing_timer, "Pack compressed");
    Utilities::pack(array, array_compressed, true);
  }
#endif

  {
    TimerOutput::Scope timer_section(computing_timer, "Pack uncompressed");
    Utilities::pack(array, array_uncompressed, false);
  }

  // computing_timer.print_summary();

  // check if compression has been invoked by comparing sizes
  deallog << "unpacked array: " << sizeof(array) << std::endl;

  deallog << "packed array without compression: " << array_uncompressed.size()
          << std::endl;

#ifdef DEAL_II_WITH_ZLIB
  deallog << "packed array with compression: " << array_compressed.size()
          << std::endl;
#endif
}


void
test()
{
  // pick large data types and arrays that could be compressed,
  // and check for both compression options
  const unsigned int N = 10000;
  double             x2[N];

  std::default_random_engine             generator(0);
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  // Test easily compressible data.
  // Default compression is much faster than best compression.
  for (unsigned int i = 0; i < N; ++i)
    {
      x2[i] = i;
    }

  check(x2);

  // Test random data, which is nearly incompressible.
  // Default compression is equally fast as best compression.
  for (unsigned int i = 0; i < N; ++i)
    {
      x2[i] = distribution(generator);
    }

  check(x2);

  deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  test();
}
