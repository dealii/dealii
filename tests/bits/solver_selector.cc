//----------------------------  solver_selector.cc  ---------------------------
//    petsc_11.cc,v 1.4 2003/07/03 10:31:46 guido Exp
//    Version: 
//
//    Copyright (C) 2004 by the deal.II authors
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

  deallog << "OK" << std::endl;
}
