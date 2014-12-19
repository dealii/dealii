// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


// check method FullMatrix::scatter_matrix_to

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>

void test ()
{
  // create a matrix with known
  // elements
  FullMatrix<double> A(5,6);
  for (unsigned int i=0; i<A.m(); ++i)
    for (unsigned int j=0; j<A.n(); ++j)
      A(i,j) = i+j;

  // pick every other row and column
  std::vector<types::global_dof_index> rows (A.m());
  for (unsigned int i=0; i<rows.size(); ++i)
    rows[i] = 2*i;

  std::vector<types::global_dof_index> cols (A.n());
  for (unsigned int i=0; i<cols.size(); ++i)
    cols[i] = 2*i;

  // do the scatter
  FullMatrix<double> X(rows.size()*2, cols.size()*2);
  A.scatter_matrix_to (rows, cols, X);

  // verify that the elements are
  // correct
  for (unsigned int i=0; i<X.m(); ++i)
    for (unsigned int j=0; j<X.n(); ++j)
      if ((i % 2 == 0) && (j % 2 == 0))
        Assert (X(i,j) == i/2 + j/2,
                ExcInternalError())
        else
          Assert (X(i,j) == 0,
                  ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main ()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);

  test();
}
