//----------------------------  conditional_ostream.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  conditional_ostream.cc  ---------------------------


// test the functions of ConditionalOStream


#include "../tests.h"
#include <base/conditional_ostream.h>
#include <fstream>
#include <iomanip>
#include <limits>


int main()
{
  std::ofstream logfile("conditional_ostream/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  ConditionalOStream o(logfile, true);
  o << "Yes" << std::endl;
  deallog << o.is_active() << std::endl;

  o.set_condition (false);
  o << "No" << std::endl;
  deallog << o.is_active() << std::endl;
}
