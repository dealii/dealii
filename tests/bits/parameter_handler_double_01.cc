// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2015 by the deal.II authors
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



// test ParameterHandler::Double description of limits

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>
#include <iomanip>


int main ()
{
  initlog();

  ParameterHandler prm;
  prm.declare_entry ("a", "1.2", Patterns::Double(), "no limit");
  prm.declare_entry ("b", "1.2", Patterns::Double(-2.13), "lower limit");
  prm.declare_entry ("c", "1.2", Patterns::Double(Patterns::Double::min_double_value, 42.0), "upper limit");
  prm.declare_entry ("d", "1.2", Patterns::Double(0.2, 42.0), "both limits");
  prm.declare_entry ("e", "1.2", Patterns::Double(1.0, -1.0), "no limits");

  prm.print_parameters (deallog.get_file_stream(), ParameterHandler::LaTeX);

  return 0;
}
