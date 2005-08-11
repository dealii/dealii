//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test functions in namespace Utilities

#include "../tests.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <base/utilities.h>


void test () 
{
  deallog << Utilities::int_to_string (42,4) << std::endl;
  deallog << Utilities::needed_digits (424) << std::endl;
  deallog << Utilities::string_to_int (" 413 ") << std::endl;

  std::vector<std::string> v;
  v.push_back ("1");
  v.push_back (" -12");
  v.push_back ("+125 ");
  Assert (Utilities::string_to_int (v).size() == 3, ExcInternalError());
  deallog << Utilities::string_to_int (v)[0] << std::endl;
  deallog << Utilities::string_to_int (v)[1] << std::endl;
  deallog << Utilities::string_to_int (v)[2] << std::endl;

  const char *p = "alpha, beta, gamma ";
  Assert (Utilities::split_comma_separated_list (p).size() == 3,
          ExcInternalError());
  Assert (Utilities::split_comma_separated_list (p)[0] == "alpha",
          ExcInternalError());  
  Assert (Utilities::split_comma_separated_list (p)[1] == "beta",
          ExcInternalError());
  Assert (Utilities::split_comma_separated_list (p)[2] == "gamma",
          ExcInternalError());

  deallog << Utilities::generate_normal_random_number (13, 44) << ' ';
  deallog << Utilities::generate_normal_random_number (13, 44) << ' ';
  deallog << Utilities::generate_normal_random_number (13, 44) << ' ';
  deallog << std::endl;
}

  
  

int main()
{
  std::ofstream logfile("utilities_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
