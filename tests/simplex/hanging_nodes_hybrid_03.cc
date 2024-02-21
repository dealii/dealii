// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// verify hanging node constraints on locally hp-refined hybrid meshes


#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/hp/fe_collection.h>

#include <vector>

#include "../tests.h"

#include "hanging_nodes.h"
#include "simplex_grids.h"


int
main(int argc, char *argv[])
{
  initlog();

  deallog.push("2d");
  {
    const unsigned int dim = 2;

    const auto cube_and_pyramid = [](Triangulation<dim> &tria) {
      GridGenerator::cube_and_pyramid(tria, 1);
    };

    for (unsigned int q = 1; q <= 4; ++q)
      for (unsigned int p = 1; p <= 2; ++p)
        {
          if (q == p)
            continue;

          deallog << "q_degree: " << q << ", p_degree: " << p << std::endl;
          test<dim>({0, 1},
                    {0, 1},
                    hp::FECollection<dim>(FE_Q<dim>(q), FE_SimplexP<dim>(p)),
                    cube_and_pyramid);
          test<dim>({1, 0},
                    {0, 1},
                    hp::FECollection<dim>(FE_Q<dim>(q), FE_SimplexP<dim>(p)),
                    cube_and_pyramid);
          test<dim>({0, 1},
                    {1, 0},
                    hp::FECollection<dim>(FE_SimplexP<dim>(p), FE_Q<dim>(q)),
                    cube_and_pyramid);
          test<dim>({1, 0},
                    {1, 0},
                    hp::FECollection<dim>(FE_SimplexP<dim>(p), FE_Q<dim>(q)),
                    cube_and_pyramid);
        }
  }
  deallog.pop();
}
