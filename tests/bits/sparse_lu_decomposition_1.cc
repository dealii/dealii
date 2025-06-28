// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// this file didn't compile at one point in time due to the private
// inheritance of SparseMatrix by SparseLUDecomposition, and the
// associated lack of accessibility of the EnableObserverPointer
// functions to the ObserverPointer
//
// it was fixed around 2003-05-22


#include <deal.II/base/observer_pointer.h>

#include <deal.II/lac/sparse_ilu.h>

#include "../tests.h"



int
main()
{
  initlog();

  ObserverPointer<SparseLUDecomposition<double>> sparse_decomp;

  deallog << "OK" << std::endl;

  return 0;
}
