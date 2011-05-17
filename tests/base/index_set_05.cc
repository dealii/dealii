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

// test IndexSet::is_contiguous and compress()

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet index_set (10);
  index_set.add_index (2);
  index_set.add_index (3);
  index_set.add_index (4);

  index_set.compress ();

  index_set.add_index (5);


  deallog << (index_set.is_contiguous() ? "true" : "false")
	  << std::endl;
  Assert (index_set.is_contiguous() == true, ExcInternalError());

  for (unsigned int i=0; i<index_set.size(); ++i)
    deallog << i << ' ' << (index_set.is_element(i) ? "true" : "false")
	    << std::endl;
}




int main()
{
  std::ofstream logfile("index_set_05/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
