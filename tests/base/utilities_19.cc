// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Utilities::inverse_Hilbert_space_filling_curve with
// effective dimension == 0

#include <deal.II/base/utilities.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;
  const std::vector<Point<dim>> points_degenerate = {Point<dim>(),
                                                     Point<dim>(),
                                                     Point<dim>()};

  const int  bit_depth = 5;
  const auto res =
    Utilities::inverse_Hilbert_space_filling_curve(points_degenerate,
                                                   bit_depth);

  for (unsigned int i = 0; i < res.size(); ++i)
    {
      for (unsigned int d = 0; d < dim; ++d)
        deallog << ' ' << res[i][d];

      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
