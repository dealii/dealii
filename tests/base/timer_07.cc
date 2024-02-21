// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/timer.h>

#include <sstream>

#include "../tests.h"

DeclException0(Timer07Exception);

int
main(int argc, char **argv)
{
  initlog();

  // capture cerr for testing purposes
  std::stringstream captured_cerr;
  std::streambuf   *old_cerr = std::cerr.rdbuf(captured_cerr.rdbuf());

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  // Test printing from TimerOutput
  try
    {
      std::cerr << "TimerOutput\n";
      TimerOutput timer_out(MPI_COMM_WORLD,
                            std::cerr,
                            TimerOutput::summary,
                            TimerOutput::cpu_times);
      timer_out.enter_subsection("Section1");

      throw Timer07Exception();
    }
  catch (const Timer07Exception &exc)
    {}

  // The Scope should still exit correctly
  try
    {
      std::cerr << "TimerOutput::Scope\n";
      TimerOutput        timer_out(MPI_COMM_WORLD,
                            std::cerr,
                            TimerOutput::summary,
                            TimerOutput::cpu_times);
      TimerOutput::Scope timer_scope(timer_out, "Section1");

      throw Timer07Exception();
    }
  catch (const Timer07Exception &exc)
    {}

  // Test that no errors are printed for MPI_COMM_SELF since no communication
  // occurs
  try
    {
      std::cerr << "TimerOutput::Scope with MPI_COMM_SELF\n";
      TimerOutput        timer_out(MPI_COMM_SELF,
                            std::cerr,
                            TimerOutput::summary,
                            TimerOutput::cpu_times);
      TimerOutput::Scope timer_scope(timer_out, "Section1");

      throw Timer07Exception();
    }
  catch (const Timer07Exception &exc)
    {}

  // convert numbers to xs to avoid printing time data
  auto        is_digit = [](const char c) -> bool { return std::isdigit(c); };
  std::string output   = captured_cerr.str();
  std::string::iterator next_number =
    std::find_if(output.begin(), output.end(), is_digit);
  while (next_number != output.end())
    {
      // convert everything between the |s to xs so that we have consistent
      // output.
      const std::string::iterator start_pipe =
        std::find(std::string::reverse_iterator(next_number),
                  output.rend(),
                  '|')
          .base();
      Assert(start_pipe != output.end(), ExcInternalError());
      const std::string::iterator end_pipe =
        std::find(next_number, output.end(), '|');
      Assert(end_pipe != output.end(), ExcInternalError());
      Assert(end_pipe - start_pipe > 1, ExcInternalError());

      std::fill(start_pipe + 1, end_pipe - 1, 'x');
      next_number = std::find_if(next_number, output.end(), is_digit);
    }


  deallog << output << std::endl;

  // restore stream
  std::cerr.rdbuf(old_cerr);
}
