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



// Test refinement of 3D wedge mesh.
// Testing for aspect ratio of refined cells

#include "deal.II/base/logstream.h"
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "./simplex_grids.h"


void
test()
{
  constexpr int      dim = 3;
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_wedges(tria, 1);

  // total number of refinements
  const unsigned int n_refinement = 4;

  for (unsigned int i = 0; i < n_refinement; i++)
    {
      FE_WedgeP<dim>   fe(1);
      MappingFE<dim>   mapping(fe);
      QGaussWedge<dim> quad(3);
      deallog << "Number of refinements: " << i << std::endl;

      // calculate the aspect ratios
      const Vector<double> aspect_ratios =
        GridTools::compute_aspect_ratio_of_cells(mapping, tria, quad);
      deallog << "aspect ratios: ";
      for (unsigned int j = 0; j < aspect_ratios.size(); j++)
        deallog << aspect_ratios[j] << " ";
      deallog << std::endl;

      // refine once for the next step
      if (i < n_refinement - 1)
        tria.refine_global(1);
    }
}

int
main()
{
  initlog();
  test();
}
