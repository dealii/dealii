// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test Support points of FE_PyramidP and FE_PyramidDGP.
#include <deal.II/fe/fe_pyramid_p.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int degree)
{
  {
    deallog << "FE_PyramidP degree: " << degree << std::endl;
    const FE_PyramidP<dim> fe(degree);
    const auto             points = fe.get_unit_support_points();

    for (const auto &p : points)
      deallog << p << std::endl;
    deallog << std::endl;

    deallog << "Unit face support points" << std::endl;
    for (const auto f : fe.reference_cell().face_indices())
      {
        deallog << "Face " << f << std::endl;
        const auto unit_face_support_points =
          fe.get_unit_face_support_points(f);
        for (const auto &p : unit_face_support_points)
          deallog << p << std::endl;
        deallog << std::endl;
      }
  }
  {
    deallog << "FE_PyramidDGP degree: " << degree << std::endl;

    const FE_PyramidDGP<dim> fe(degree);
    const auto               points = fe.get_unit_support_points();

    for (const auto &p : points)
      deallog << p << std::endl;
    deallog << std::endl;
  }
}

int
main()
{
  initlog();

  for (unsigned int i = 1; i < 4; ++i)
    test<3>(i);
}
