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


// tests for the ComponentMask class
//
// here: test that creating a mask with constant elements using the direct
// constructor


#include "../tests.h"
#include <deal.II/fe/component_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  std::vector<bool> v(12, false);
  ComponentMask     m(12, false);

				   // verify equality
  for (unsigned int i=0; i<v.size(); ++i)
    Assert (m[i] == v[i], ExcInternalError());

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile ("component_mask_03/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
