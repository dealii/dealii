// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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
