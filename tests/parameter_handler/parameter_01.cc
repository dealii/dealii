// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// check the ParameterHandler::add_parameter() function


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>


template <typename T>
void check (const T            initializer,
            const std::string &other_value)
{
  T parameter = initializer;

  ParameterHandler prm;
  prm.add_parameter ("parameter", parameter, "", Patterns::DefaultPattern<T>::create());

  // see what the default value of the parameter is
  deallog << "default=" << prm.get ("parameter") << std::endl;

  // now read a different value and output the value so read again
  const std::string p = "set parameter = " + other_value + "\n";
  std::istringstream in(p);
  prm.parse_input (in);

  deallog << "new value=" << prm.get ("parameter") << std::endl;
}



int main ()
{
  initlog();

  check<int> (13, "42");
  check<std::vector<std::string>> ({ "1","2","3" }, "4,5,6");
}
