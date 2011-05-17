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

// test IndexSet::operator &

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet is1 (100);
  IndexSet is2 (100);

				   // randomly add 90 elements to each
				   // set, some of which may be
				   // repetitions of previous ones
  for (unsigned int i=0; i<9*is1.size()/10; ++i)
    {
      is1.add_index (rand() % is1.size());
      is2.add_index (rand() % is2.size());
    }

  IndexSet is3 = is1 & is2;

  deallog << "Set sizes: "
	  << is1.n_elements() << ' '
    	  << is2.n_elements() << ' '
    	  << is3.n_elements() << std::endl;

  for (unsigned int i=0; i<is3.size(); ++i)
    {
      deallog << i << ' ' << (is3.is_element(i) ? "true" : "false")
	      << std::endl;

      Assert ((is1.is_element(i) && is2.is_element(i))
	      ==
	      is3.is_element(i),
	      ExcInternalError());
    }

				   // some sanity tests
  Assert ((is1 & is2) == (is2 & is1), ExcInternalError());
  Assert ((is1 & is3) == (is2 & is3), ExcInternalError());
  Assert ((is1 & is3) == is3, ExcInternalError());
  Assert ((is3 & is1) == is3, ExcInternalError());
  Assert ((is3 & is3) == is3, ExcInternalError());

  deallog << "OK" << std::endl;
}




int main()
{
  std::ofstream logfile("index_set_14/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
