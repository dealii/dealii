// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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
