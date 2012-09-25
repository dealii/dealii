//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// tests for the BlockMask class
//
// here: test that creating a mask from a vector<bool> works


#include "../tests.h"
#include <deal.II/fe/block_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  std::vector<bool> v(12);
  for (unsigned int i=0; i<v.size(); ++i)
    v[i] = (i%3 == 0);

  BlockMask m(v);

				   // verify equality
  for (unsigned int i=0; i<v.size(); ++i)
    Assert (m[i] == v[i], ExcInternalError());

				   // this needs to throw an exception
  m[v.size()];

  deallog << "OK" << std::endl;
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();
  std::ofstream logfile ("block_mask_02/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
