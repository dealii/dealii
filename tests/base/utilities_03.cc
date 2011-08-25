//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test functions in namespace Utilities

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/utilities.h>


void test ()
{
  deallog << Utilities::string_to_double (" 413 ") << std::endl;

  std::vector<std::string> v;
  v.push_back ("1.5");
  v.push_back (" -12.5");
  v.push_back ("+125.5 ");
  Assert (Utilities::string_to_double (v).size() == 3, ExcInternalError());
  deallog << Utilities::string_to_double (v)[0] << std::endl;
  deallog << Utilities::string_to_double (v)[1] << std::endl;
  deallog << Utilities::string_to_double (v)[2] << std::endl;
}




int main()
{
  std::ofstream logfile("utilities_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
