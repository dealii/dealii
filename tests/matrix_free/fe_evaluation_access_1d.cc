// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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
