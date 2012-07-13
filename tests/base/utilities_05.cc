//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// Utilities::get_integer_at_position

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>

#include <deal.II/base/utilities.h>

using namespace dealii;




void test ()
{
  int number = 5;
  for (unsigned int i=0; i<7; ++i)
    {
      std::ostringstream s;
      s << "test test" << number << "test test";

      Assert (Utilities::get_integer_at_position (s.str(),
						  9).first
	      == number,
	      ExcInternalError());
      Assert (Utilities::get_integer_at_position (s.str(),
						  9).second
	      == i+1,
	      ExcInternalError());

      deallog << i << ' ' << Utilities::get_integer_at_position (s.str(),
								 9).first
	      << std::endl;

      number = number*10 + i;
    }
}




int main()
{
  std::ofstream logfile("utilities_05/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
