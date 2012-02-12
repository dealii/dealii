//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test IndexSet::fill_binary_vector

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

				   // randomly add 90 elements to each
				   // set, some of which may be
				   // repetitions of previous ones
  for (unsigned int i=0; i<9*is1.size()/10; ++i)
    is1.add_index (rand() % is1.size());

  std::vector<bool> zeros_and_ones(is1.size());
  is1.fill_binary_vector (zeros_and_ones);

  deallog << "Original index set: " << std::endl;
  is1.print(deallog);

  for (unsigned int i=0; i<is1.size(); i++)
    Assert(is1.is_element(i) == zeros_and_ones[i],
	   ExcInternalError());

  deallog << "OK" << std::endl;
}




int main()
{
  std::ofstream logfile("index_set_22/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
