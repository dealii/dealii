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
