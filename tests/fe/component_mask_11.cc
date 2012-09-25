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
// here: ComponentMask::operator==


#include "../tests.h"
#include <deal.II/fe/component_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  std::vector<bool> v1(12);
  for (unsigned int i=0; i<v1.size(); ++i)
    v1[i] = (i%3 == 0);
  std::vector<bool> v2(12);
  for (unsigned int i=0; i<v2.size(); ++i)
    v2[i] = (i%4 == 0);

  std::vector<bool> v(12);
  for (unsigned int i=0; i<v.size(); ++i)
    v[i] = (v1[i] || v2[i]);

  ComponentMask m1(v1);
  ComponentMask m2(v2);
  ComponentMask m = m1 | m2;

				   // verify equality
  Assert (m == ComponentMask(v),
	  ExcInternalError());
  Assert (!(m == m1),
	  ExcInternalError());
  Assert (!(m == ComponentMask(v1)),
	  ExcInternalError());
  Assert (!(m == ComponentMask(v2)),
	  ExcInternalError());

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile ("component_mask_11/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
