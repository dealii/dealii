// ---------------------------------------------------------------------
//
// Copyright (C) 2016-2022 by the deal.II authors
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
namespace callgrind_wrapper
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
      const auto my_pid = getpid(); // do something useful
      CALLGRIND_DUMP_STATS_AT("callgrind-wrapper-token");
      CALLGRIND_STOP_INSTRUMENTATION;

      std::ifstream      callgrind_output("callgrind.out");
      std::string        token;
      unsigned long long cycles      = 0ull;
      bool               found_token = false;
      unsigned long      found_pid   = 0;

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
        "callgrind_wrapper::start_instrumentation() can only be called when "
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
  inline unsigned long long
  stop_instrumentation()
  {
    CALLGRIND_DUMP_STATS;
    CALLGRIND_STOP_INSTRUMENTATION;

    std::ifstream      callgrind_output("callgrind.out");
    std::string        token;
    unsigned long long cycles = 0ull;

    while (callgrind_output)
      {
        callgrind_output >> token;
        if (token == "summary:")
          callgrind_output >> cycles;
      }

    std::remove("callgrind.out");
    return cycles;
  }

} // namespace callgrind_wrapper

#endif
