// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test refinement of 3D simplex mesh.
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



void
test()
{
  constexpr int      dim = 3;
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1, 0, 1);

  for (unsigned int i = 0; i < 3; i++)
    {
      FE_SimplexP<dim>   fe(1);
      MappingFE<dim>     mapping(fe);
      QGaussSimplex<dim> quad(3);

      const Vector<double> aspect_ratios =
        GridTools::compute_aspect_ratio_of_cells(mapping, tria, quad);
      deallog << "Number of refinements: " << i << std::endl;
      deallog << "aspect ratios: ";
      for (unsigned int j = 0; j < aspect_ratios.size(); j++)
        deallog << aspect_ratios[j] << " ";
      deallog << std::endl;
      if (i < 2)
        tria.refine_global(1);
    }
}

int
main()
{
  initlog();
  test();
}
