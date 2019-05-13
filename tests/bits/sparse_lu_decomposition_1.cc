// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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


// this file didn't compile at one point in time due to the private
// inheritance of SparseMatrix by SparseLUDecomposition, and the
// associated lack of accessibility of the Subscriptor functions to
// the SmartPointer
//
// it was fixed around 2003-05-22


#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/sparse_ilu.h>

#include "../tests.h"



int
main()
{
  initlog();

  SmartPointer<SparseLUDecomposition<double>> sparse_decomp;

  deallog << "OK" << std::endl;

  return 0;
}
