//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test IndexSet::clear

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet index_set (20);
  index_set.add_range (2,4);
  index_set.add_range (12,18);
  index_set.add_index (6);
  index_set.add_index (8);
  index_set.add_index (14);
  index_set.add_index (16);

				   // clear the IndexSet and then set elements
				   // again
  index_set.clear ();

  index_set.add_range (2,4);
  index_set.add_range (12,18);
  index_set.add_index (6);
  index_set.add_index (8);
  index_set.add_index (14);
  index_set.add_index (16);

  for (unsigned int i=0; i<index_set.size(); ++i)
    deallog << i << ' ' << (index_set.is_element(i) ? "true" : "false")
	    << std::endl;
}




int main()
{
  std::ofstream logfile("index_set_21/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
