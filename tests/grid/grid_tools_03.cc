// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check GridTools::volume


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test1()
{
  // test 1: hypercube
  if (true)
    {
      Triangulation<dim> tria;
      GridGenerator::hyper_cube(tria);

      for (unsigned int i = 0; i < 2; ++i)
        {
          tria.refine_global(2);
          deallog << dim << "d, "
                  << "hypercube volume, " << i * 2
                  << " refinements: " << GridTools::volume(tria) << std::endl;
        }

      Triangulation<dim> simplex_tria;
      GridGenerator::convert_hypercube_to_simplex_mesh(tria, simplex_tria);
      deallog << dim << "d, simplex hypercube volume: "
              << GridTools::volume(simplex_tria) << std::endl;
    }

  // test 2: hyperball
  if (dim >= 2)
    {
      Triangulation<dim> tria;
      GridGenerator::hyper_ball(tria, Point<dim>(), 1);

      for (unsigned int i = 0; i < 4; ++i)
        {
          tria.refine_global(1);
          deallog << dim << "d, "
                  << "hyperball volume, " << i
                  << " refinements: " << GridTools::volume(tria) << std::endl;
        }

      Triangulation<dim> simplex_tria;
      GridGenerator::convert_hypercube_to_simplex_mesh(tria, simplex_tria);
      simplex_tria.set_all_manifold_ids(numbers::flat_manifold_id);
      deallog << dim << "d, simplex hyperball volume: "
              << GridTools::volume(simplex_tria) << std::endl;

      deallog << "exact value="
              << (dim == 2 ? numbers::PI : 4. / 3. * numbers::PI) << std::endl;
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test1<1>();
  test1<2>();
  test1<3>();

  return 0;
}
