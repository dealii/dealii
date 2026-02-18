// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check that VectorTools::project works for FESystem(FE_Q) elements correctly
// on a uniformly refined mesh for functions of degree q

#include "project_common.h"


template <int dim>
void
test()
{
  for (unsigned int p = 1; p < 5 - dim; ++p)
    test_with_wrong_face_orientation(
      FESystem<dim>(FE_Q<dim>(p), 1, FE_DGQ<dim>(p + 1), 1), p);
}
