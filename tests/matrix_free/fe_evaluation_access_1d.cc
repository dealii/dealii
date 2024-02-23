// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// FEEvaluationAccess<1,1,double> didn't compile because we had
// conflicting partial specializations of this class

#include <deal.II/matrix_free/fe_evaluation.h>

#include "../tests.h"


int
main()
{
  initlog();

  FEEvaluationAccess<1, 1, double, false> *test; // didn't compile before
  deallog << "OK" << std::endl;
}
