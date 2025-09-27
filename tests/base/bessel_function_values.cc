// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify that the new C++17 Bessel function produces equivalent output to jn.


#include <deal.II/base/function_bessel.h>

#include "../tests.h"


int
main()
{
  initlog();

  for (unsigned int order = 0; order < 3; order += 1)
    for (double wave_n = 0.; wave_n < 4.0; wave_n += 1.5)
      for (double x = 1.0; x < 2.0; x += 0.5)
        for (double y = 2.0; y < 3.0; y += 0.5)
          {
            Functions::Bessel1<2> bessel(order, wave_n, Point<2>(x, y));

            const Point<2> new_point(x + 0.1, y + 0.1);
            deallog << "value at (" << new_point[0] << ", " << new_point[1]
                    << ") = " << bessel.value(new_point) << std::endl;
            deallog << "gradient at (" << new_point[0] << ", " << new_point[1]
                    << ") = " << bessel.gradient(new_point) << std::endl;
          }
}
