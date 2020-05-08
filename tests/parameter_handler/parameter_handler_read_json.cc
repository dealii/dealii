// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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



// check ParameterHandler::parse_input_from_xml

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  // default values
  int         int1         = 1;
  int         int2         = 2;
  double      double1      = 1.234;
  double      double2      = 4.321;
  std::string str          = "< & > ; /";
  int         intint       = 2;
  double      doubledouble = 6.1415926;

  ParameterHandler prm;
  prm.add_parameter("int1", int1, "doc 1");
  prm.add_parameter("int2", int2, "doc 2");
  prm.enter_subsection("ss1");
  {
    prm.add_parameter("double 1", double1, "doc 3");

    prm.enter_subsection("ss2");
    {
      prm.add_parameter("double 2", double2, "doc 4");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();

  // things with strange characters
  prm.enter_subsection("Testing%testing");
  {
    prm.add_parameter("string&list", str, "docs 1");
    prm.add_parameter("int*int", intint);
    prm.add_parameter("double+double", doubledouble, "docs 3");
  }
  prm.leave_subsection();

  // read from json
  std::ifstream in(SOURCE_DIR "/prm/parameter_handler_read_json.prm");
  prm.parse_input_from_json(in);

  Assert(int1 == 2, ExcNotImplemented());
  Assert(int2 == 3, ExcNotImplemented());
  Assert(double1 == 2.234, ExcNotImplemented());
  Assert(double2 == 5.321, ExcNotImplemented());
  Assert(str == "< & > ; /", ExcNotImplemented());
  Assert(intint == 2, ExcNotImplemented());
  Assert(doubledouble == 7.1415926, ExcNotImplemented());

  // write it out again
  prm.print_parameters(deallog.get_file_stream(), ParameterHandler::JSON);
  deallog.get_file_stream() << std::endl;

  return 0;
}
