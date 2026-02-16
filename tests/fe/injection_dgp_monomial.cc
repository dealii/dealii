// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



char logname[] = "output";


#include "injection_common.h"


template <int dim>
void
test()
{
  for (unsigned int i = 1; i < 4; ++i)
    for (unsigned int j = i; j < 4; ++j)
      do_check(FE_DGPMonomial<dim>(i), FE_DGPMonomial<dim>(j));
}
