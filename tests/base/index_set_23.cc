//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// we managed to get a function as simple as IndexSet::operator== wrong -- who
// knew?

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet is1 (100), is2 (100);

  is1.add_range (0,10);
  is2.add_range (0,20);

  Assert((is1 == is2) == false, ExcInternalError());
  Assert((is1 != is2) == true, ExcInternalError());

  deallog << "OK" << std::endl;
}




int main()
{
  std::ofstream logfile("index_set_23/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
