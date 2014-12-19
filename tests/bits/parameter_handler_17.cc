// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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



// One can (accidentally) call ParameterHandler::get_int/double on parameters
// that are really strings. There used to be a bug in these functions in that
// they didn't throw an error when the string wasn't actually convertible to a
// number

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>

void check ()
{
  ParameterHandler prm;
  prm.declare_entry ("a", "this that and the other", Patterns::Anything(),
		     "");
  try
    {
      prm.get_double ("a");
    }
  catch (const ParameterHandler::ExcConversionError &)
    {
      deallog << "get_double() detected the mistake" << std::endl;
    }

  try
    {
      prm.get_integer ("a");
    }
  catch (const ParameterHandler::ExcConversionError &)
    {
      deallog << "get_integer() detected the mistake" << std::endl;
    }

  try
    {
      prm.get_bool ("a");
    }
  catch (const ParameterHandler::ExcConversionError &)
    {
      deallog << "get_bool() detected the mistake" << std::endl;
    }
  
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check ();

  return 0;
}
