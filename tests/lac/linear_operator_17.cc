// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test the LinearOperator template on a trivial vector implementation
// :: RightVector -> LeftVector

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"



int
main()
{
  initlog();

  const std::function<void(Vector<double> &, bool)> reinit_vector =
    [](Vector<double> &v, bool omit_zeroing_entries) {
      v.reinit(3, omit_zeroing_entries);
    };

  const auto filter = mean_value_filter(reinit_vector);

  Vector<double> vec(3);
  vec[0] = 1.;
  vec[1] = 2.;
  vec[2] = 3.;

  deallog << vec << std::endl;

  filter.vmult_add(vec, vec);
  deallog << vec << std::endl;

  filter.vmult(vec, vec);
  deallog << vec << std::endl;
}
