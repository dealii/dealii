//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test IndexSet::n_elements()

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <base/index_set.h>


void test ()
{
  IndexSet index_set (20);
  index_set.add_index (2);
  index_set.add_index (3);
  index_set.add_index (4);

  index_set.add_index (6);
  index_set.add_index (7);

  index_set.add_index (9);

  deallog << index_set.n_elements()
	  << std::endl;
}




int main()
{
  std::ofstream logfile("index_set_09/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
