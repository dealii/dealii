// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// SparsityPattern::copy_from crashed when the number of rows or columns
// was zero

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <list>
#include <set>

#include "../testmatrix.h"
#include "../tests.h"


int
main()
{
  std::ofstream logfile("output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  DynamicSparsityPattern csp(10, 0);

  SparsityPattern sp;
  sp.copy_from(csp);

  deallog << "OK" << std::endl;
}
