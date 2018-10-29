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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test that the wall time calculated by a timer agrees with the value
// calculated by the system clock. This verifies that we fixed a bug where the
// wall time was doubled.

#include <deal.II/base/timer.h>

#include <chrono>
#include <thread>

#include "../tests.h"

int
main(int argc, char **argv)
{
  initlog();

  Timer timer;
  timer.start();
  const auto t0 = std::chrono::system_clock::now();
  std::this_thread::sleep_for(std::chrono::seconds(4));
  timer.stop();
  const auto t1 = std::chrono::system_clock::now();

  // verify that the timer wall time is not double the manually calculated one
  AssertThrow(std::abs(
                double(std::chrono::duration_cast<std::chrono::seconds>(t1 - t0)
                         .count()) -
                timer.wall_time()) < 0.5,
              ExcMessage("The measured times should be close."));

  deallog << "OK" << std::endl;
}
