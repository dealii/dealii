// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// Oliver found an example, where sparse_matrix_iterator->value=0 didn't work,
// because the iterator->value expects a double on the right hand side, not an
// integer. If the right hand side is zero, it can also be converted to a
// pointer, which leads to an ambiguity. Fix this by having an additional
// operator= in the iterator/reference class

#include "../tests.h"
#include <deal.II/lac/sparse_matrix.h>
#include <iomanip>
#include <fstream>


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // this test only needs to compile, not run
  if (false)
    {
      SparseMatrix<double>::iterator *i;
      (*i)->value () = (int)0;
    }

  deallog << "OK" << std::endl;

  return 0;
}
