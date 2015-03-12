// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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


// test Utilities::split_string_list with a string that only contains
// the delimiter and, possibly, spaces

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/utilities.h>


void test ()
{
  // verify the documented behavior of eating trailing delimiters
  {
    deallog << Utilities::split_string_list (",").size()
	    << std::endl;
    deallog << Utilities::split_string_list (" , ").size()
	    << std::endl;
  }

  {
    deallog << Utilities::split_string_list (",,").size()
	    << std::endl;
    deallog << Utilities::split_string_list (" , , ").size()
	    << std::endl;
  }

  // try some more esoteric cases:
  {
    deallog << Utilities::split_string_list (" , , ", ' ').size()
	    << std::endl;
  }

  {
    deallog << Utilities::split_string_list (" ", ' ').size()
	    << std::endl;
    deallog << Utilities::split_string_list ("   ", ' ').size()
	    << std::endl;
  }

  Assert (Utilities::split_string_list(" ; ", ';').size() == 1,
	  ExcInternalError());
  Assert (Utilities::split_string_list(" ; ", ';')[0] == "",
	  ExcInternalError());
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
