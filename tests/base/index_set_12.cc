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

// test indexing in IndexSet variables for a non-contiguous range

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet index_set (20);

  index_set.add_index (2);
  index_set.add_index (3);
  index_set.add_index (4);

  index_set.add_index (6);
  index_set.add_index (7);

  index_set.compress ();

  index_set.add_index (9);

  for (unsigned int i=0; i<index_set.n_elements(); ++i)
    {
      deallog << index_set.nth_index_in_set(i)
	      << std::endl;
      Assert (index_set.index_within_set(index_set.nth_index_in_set(i)) == i,
	      ExcInternalError());
    }
  deallog << "OK" << std::endl;

  for (unsigned int i=0; i<index_set.size(); ++i)
    if (index_set.is_element (i))
      deallog << i << ' ' << index_set.index_within_set(i)
	      << std::endl;
}




int main()
{
  std::ofstream logfile("index_set_12/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
