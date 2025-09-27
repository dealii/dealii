// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



char logname[] = "output";

#include <deal.II/base/quadrature_lib.h>

#include "injection_common.h"


template <int dim>
void
test()
{
  deallog << std::setprecision(8);

  for (unsigned int i = 1; i < 4; ++i)
    for (unsigned int j = i; j < 4; ++j)
      do_check(FESystem<dim>(FE_Q<dim>(QIterated<1>(QTrapezoid<1>(), i)),
                             1,
                             FE_DGQ<dim>(i - 1),
                             1),
               FESystem<dim>(FE_Q<dim>(QIterated<1>(QTrapezoid<1>(), j)),
                             1,
                             FE_DGQ<dim>(j - 1),
                             1));
}
