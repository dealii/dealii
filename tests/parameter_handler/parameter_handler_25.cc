// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <map>

#include "../tests.h"

void
success()
{
  unsigned int dim       = 2;
  std::string  precision = "double";

  ParameterHandler prm;
  prm.enter_subsection("General");
  // this one does not have to be set
  prm.add_parameter("dim",
                    dim,
                    "Number of space dimensions",
                    Patterns::Integer(2, 3));
  // this one has to be set
  prm.declare_entry("Precision",
                    precision,
                    Patterns::Selection("float|double"),
                    "Floating point precision",
                    true);
  prm.leave_subsection();

  // this try-catch simulates parsing from an incomplete/incorrect input file
  try
    {
      prm.enter_subsection("General");
      prm.set("Precision", "float");
      prm.leave_subsection();
    }
  catch (const std::exception &exc)
    {
      deallog << exc.what() << std::endl;
    }

  // check set status
  try
    {
      prm.assert_that_entries_have_been_set();
    }
  catch (const std::exception &exc)
    {
      deallog << exc.what() << std::endl;
    }

  deallog << std::endl << "successful" << std::endl;
}

void
fail()
{
  unsigned int dim       = 2;
  std::string  precision = "double";

  ParameterHandler prm;
  prm.enter_subsection("General");
  // both parameters have to be set
  prm.add_parameter(
    "dim", dim, "Number of space dimensions", Patterns::Integer(2, 3), true);
  prm.add_parameter("Precision",
                    precision,
                    "Floating point precision",
                    Patterns::Selection("float|double"),
                    true);
  prm.leave_subsection();

  // this try-catch simulates parsing from an incomplete/incorrect input file
  try
    {
      prm.enter_subsection("General");
      prm.set("Precison", "float"); // here is a typo!
      // dim is not set!
      prm.leave_subsection();
    }
  catch (const std::exception &exc)
    {
      // Print everything after the "violated condition" part
      // of the error message, assuming the condition is shown.
      // If it isn't (because the condition was simply 'false' and
      // so the error printing part suppresses this part), then
      // show the part after "Additional information":
      std::string error = exc.what();
      if (auto start = error.find("The violated condition was:");
          start != std::string::npos)
        deallog << error.substr(start) << std::endl;
      else if (auto start = error.find("Additional information:");
               start != std::string::npos)
        deallog << error.substr(start) << std::endl;
    }

  // check set status
  try
    {
      prm.assert_that_entries_have_been_set();
    }
  catch (const std::exception &exc)
    {
      // Same as above:
      std::string error = exc.what();
      if (auto start = error.find("The violated condition was:");
          start != std::string::npos)
        deallog << error.substr(start) << std::endl;
      else if (auto start = error.find("Additional information:");
               start != std::string::npos)
        deallog << error.substr(start) << std::endl;
    }
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  try
    {
      success();
      fail();
    }
  catch (const std::exception &exc)
    {
      deallog << exc.what() << std::endl;
    }
}
