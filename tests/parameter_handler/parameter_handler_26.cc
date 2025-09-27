// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
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
success(const std::string &filename)
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

  try
    {
      prm.parse_input(filename, "", true, true);
    }
  catch (const std::exception &exc)
    {
      deallog << exc.what() << std::endl;
    }

  deallog << std::endl << "successful" << std::endl;
}

void
fail(const std::string &filename)
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

  try
    {
      prm.parse_input(filename, "", true, true);
    }
  catch (const std::exception &exc)
    {
      deallog << exc.what() << std::endl;
    }
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  const std::string source = SOURCE_DIR;
  try
    {
      success(source + "/prm/parameter_handler_26_success.json");
      fail(source + "/prm/parameter_handler_26_fail.json");
      success(source + "/prm/parameter_handler_26_success.prm");
      fail(source + "/prm/parameter_handler_26_fail.prm");
    }
  catch (const std::exception &exc)
    {
      deallog << exc.what() << std::endl;
    }
}
