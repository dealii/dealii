//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test Utilities::fixed_power

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/utilities.h>


template <int dim>
void test () 
{
  deallog << Utilities::fixed_power<dim> (2) << std::endl;
  deallog << Utilities::fixed_power<dim> (-2) << std::endl;
  deallog << Utilities::fixed_power<dim> (2.5) << std::endl;
  deallog << Utilities::fixed_power<dim> (-2.5) << std::endl;
  deallog << std::endl;
}

  
  

int main()
{
  std::ofstream logfile("utilities_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  test<4> ();
}
