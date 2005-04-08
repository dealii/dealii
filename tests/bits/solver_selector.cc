//----------------------------  solver_selector.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver_selector.cc  ---------------------------

// I managed to break solver selector once by making some variables in the
// AdditionalData structures of the solvers const. This test simply
// instantiates that class, to make sure it still compiles

#include <lac/solver_selector.h>
#include <fstream>

// instantiation here
template class SolverSelector<>;

int main () 
{  
  std::ofstream logfile("solver_selector.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << "OK" << std::endl;
}
