// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

  try
    {
      std::string source   = SOURCE_DIR;
      std::string filename = source + "/prm/parameter_handler_26_success.json";
      prm.parse_input(filename, "", true, true);
    }
  catch (std::exception &exc)
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

  try
    {
      std::string source   = SOURCE_DIR;
      std::string filename = source + "/prm/parameter_handler_26_fail.json";
      prm.parse_input(filename, "", true, true);
    }
  catch (std::exception &exc)
    {
      deallog << exc.what() << std::endl;
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
  catch (std::exception &exc)
    {
      deallog << exc.what() << std::endl;
    }
}
