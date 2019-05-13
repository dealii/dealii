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



//   Like parameter_handler_03, but use a MultipleSelection pattern that starts
//   with a space; eat that space

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  try
    {
      initlog();

      ParameterHandler prm;
      prm.enter_subsection("Testing");
      prm.declare_entry("string list1",
                        "a",
                        Patterns::List(Patterns::Selection(" a|b|c|d|e|f|g|h")),
                        "docs 1");
      prm.declare_entry("string list2",
                        "h",
                        Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h ")),
                        "docs 2");
      prm.declare_entry("int", "1", Patterns::Integer());
      prm.declare_entry("double", "3.1415926", Patterns::Double(), "docs 3");
      prm.leave_subsection();

      // read and then write parameters
      prm.parse_input(SOURCE_DIR "/prm/parameter_handler_3_with_space.prm");
      prm.print_parameters(deallog.get_file_stream(), ParameterHandler::Text);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };

  return 0;
}
