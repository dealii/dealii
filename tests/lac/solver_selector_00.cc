// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// I managed to break solver selector once by making some variables in the
// AdditionalData structures of the solvers const. This test simply
// instantiates that class, to make sure it still compiles

#include "../tests.h"
#include <deal.II/lac/solver_selector.h>
#include <fstream>

DEAL_II_NAMESPACE_OPEN
// instantiation here
template class SolverSelector<>;
DEAL_II_NAMESPACE_CLOSE

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << "OK" << std::endl;
}
