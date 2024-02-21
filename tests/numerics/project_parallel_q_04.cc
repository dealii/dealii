// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that VectorTools::project works for Q elements correctly

#include "project_parallel_common.h"


template <int dim>
void
test()
{
  test_with_2d_deformed_mesh<dim, 2, 1>(FESystem<dim>(FE_Q<dim>(1), 2));
  test_with_2d_deformed_mesh<dim, 1, 2>(FESystem<dim>(FE_Q<dim>(2), 1));
}
