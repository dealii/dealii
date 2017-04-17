// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check the ParameterHandler::add_action() function. like _01, but
// attach a second action


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>


void check (const char *p)
{
  std::string parameter_set_by_action;

  ParameterHandler prm;
  prm.declare_entry ("test_1", "-1,0",
                     Patterns::List(Patterns::Integer(-1,1),2,3));
  prm.add_action ("test_1",
                  [&](const std::string &s)
  {
    deallog << "In action 1:" << s << std::endl;
    parameter_set_by_action = s;
    return true;
  });
  prm.add_action ("test_1",
                  [&](const std::string &s)
  {
    deallog << "In action 2:" << s << std::endl;
    parameter_set_by_action = s + " some modification";
    return true;
  });


  std::ifstream in(p);

  deallog << "Reading parameters" << std::endl;
  prm.read_input (in);

  deallog << "test_1=" << prm.get ("test_1") << std::endl;
  deallog << "Saved parameter: " << parameter_set_by_action << std::endl;
}



int main ()
{
  initlog();

  check (SOURCE_DIR "/prm/parameter_handler_1.prm");

  return 0;
}
