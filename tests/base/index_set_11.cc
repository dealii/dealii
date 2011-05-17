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

// test IndexSet::add_indices

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet index_set (10);

  {
    const unsigned int array[] = { 2, 3, 4 };
    index_set.add_indices ((const unsigned int*)array, array+sizeof(array)/sizeof(array[0]));
  }
  {
    const unsigned int array[] = { 7,6 };
    index_set.add_indices ((const unsigned int*)array, array+sizeof(array)/sizeof(array[0]));
  }

  Assert (index_set.is_contiguous() == false, ExcInternalError());

  for (unsigned int i=0; i<index_set.size(); ++i)
    deallog << i << ' ' << (index_set.is_element(i) ? "true" : "false")
	    << std::endl;


  {
    const unsigned int array[] = { 5 };
    index_set.add_indices ((const unsigned int*)array, array+sizeof(array)/sizeof(array[0]));
  }

  Assert (index_set.is_contiguous() == true, ExcInternalError());

  for (unsigned int i=0; i<index_set.size(); ++i)
    deallog << i << ' ' << (index_set.is_element(i) ? "true" : "false")
	    << std::endl;


  deallog << "OK" << std::endl;
}




int main()
{
  std::ofstream logfile("index_set_11/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
