// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022-2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_valgrind_instrumentation_h
#define dealii_valgrind_instrumentation_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <unistd.h>
#include <valgrind/callgrind.h>

#include <cstdio>
#include <fstream>
#include <mutex>


/**
 * A simple wrapper around callgrind to query current instruction counter
 * values. In order to make use of the functions in this namespace the
 * executable has to be started with valgrind:
 * @code
 *   valgrind --tool=callgrind -q --combine-dumps=yes --instr-atstart=no \
 *            --callgrind-out-file=callgrind.out ./executable
 * @endcode
 */
namespace CallgrindWrapper
{
  /**
   * Reset the current instruction counter to zero and start callgrind
   * instrumentation.
   */
  DEAL_II_ALWAYS_INLINE
  inline void
  start_instrumentation()
  {
    /*
     * Verify that we are running under valgrind and that valgrind is
     * called with the right parameters.
     */
    static bool is_initialized = []() {
      // Is this executable running under valgrind supervision?
      if (RUNNING_ON_VALGRIND != 1)
        return false;

      // Start and stop valgrind instrumentation
      CALLGRIND_STOP_INSTRUMENTATION;
      std::remove("callgrind.out");
      CALLGRIND_ZERO_STATS;
      CALLGRIND_START_INSTRUMENTATION;
      const unsigned my_pid = getpid(); // do something useful
      CALLGRIND_DUMP_STATS_AT("callgrind-wrapper-token");
      CALLGRIND_STOP_INSTRUMENTATION;

      std::ifstream callgrind_output("callgrind.out");
      std::string   token;
      std::uint64_t cycles      = 0ull;
      bool          found_token = false;
      unsigned long found_pid   = 0;

      while (callgrind_output)
        {
          callgrind_output >> token;
          if (token == "pid:")
            callgrind_output >> found_pid;
          if (token == "desc:")
            {
              std::string line;
              std::getline(callgrind_output, line);
              if (line == " Trigger: Client Request: callgrind-wrapper-token")
                found_token = true;
            }
          if (token == "summary:")
            callgrind_output >> cycles;
        }

      if (!found_token)
        return false;

      if (found_pid != my_pid)
        return false;

      if (cycles == 0)
        return false;

      return true;
    }();

    AssertThrow(
      is_initialized,
      dealii::ExcMessage(
        "CallgrindWrapper::start_instrumentation() can only be called when "
        "the executable is run via \"valgrind --tool=callgrind -q "
        "--combine-dumps=yes --instr-atstart=no "
        "--callgrind-out-file=callgrind.out ./executable\""));

    CALLGRIND_ZERO_STATS;
    CALLGRIND_START_INSTRUMENTATION;
  }


  /**
   * Stop callgrind instrumentation and return the number of instructions
   * executed since the last start_instrumentation() call.
   */
  DEAL_II_ALWAYS_INLINE
  inline std::uint64_t
  stop_instrumentation()
  {
    CALLGRIND_DUMP_STATS;
    CALLGRIND_STOP_INSTRUMENTATION;

    std::ifstream callgrind_output("callgrind.out");
    std::string   token;
    std::uint64_t cycles = 0ull;

    while (callgrind_output)
      {
        callgrind_output >> token;
        if (token == "summary:")
          callgrind_output >> cycles;
      }

    std::remove("callgrind.out");
    return cycles;
  }



  /**
   * A function that counts the cycles necessary to execute the given function
   * argument. The following are therefore equivalent:
   * @code
   *   start_instrumentation();
   *   my_function();
   *   const auto cycles = stop_instrumentation();
   * @endcode
   * and
   * @code
   *   const auto cycles = count_cycles([](){ my_function(); });
   * @endcode
   * If the call to `my_function()` involves arguments or code pieces that
   * access variables in the environment of the place where the function is
   * called, then one can of course capture these variables in the construction
   * of the lambda function and instead call variations such as
   * @code
   *   const auto cycles = count_cycles([&](){ my_function(); });
   * @endcode
   */
  template <typename Func>
  inline std::uint64_t
  count_cycles(Func &&f)
  {
    start_instrumentation();
    f();
    return stop_instrumentation();
  }
} // namespace CallgrindWrapper

#endif
