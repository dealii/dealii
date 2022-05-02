// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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



// verify hanging node constraints on locally h-refined hybrid meshes
/*
 * dofs will be enumerated as follows with d = 1
 *  scenario 1:        scenario 2:        scenario 3:
 *   2-------3          2-------3          7---8---9
 *   |       |\         |       |\         |   |   |\
 *   |       | \        |       | 6        |   |   | \
 *   |       |  \       |       |/|\       |   |   |  \
 *   |       |   4      |       4 | 7      3---4---6   0
 *   |       |  /       |       |\|/       |   |   |  /
 *   |       | /        |       | 5        |   |   | /
 *   |       |/         |       |/         |   |   |/
 *   0-------1          0-------1          1---2---5
 */


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

    for (unsigned int d = 1; d <= 2; ++d)
      {
        deallog << "degree: " << d << std::endl;
        test<dim>({0, 1},
                  {0, 1},
                  hp::FECollection<dim>(FE_Q<dim>(d), FE_SimplexP<dim>(d)),
                  cube_and_pyramid);
        test<dim>({1, 0},
                  {0, 1},
                  hp::FECollection<dim>(FE_Q<dim>(d), FE_SimplexP<dim>(d)),
                  cube_and_pyramid);
        test<dim>({0, 1},
                  {1, 0},
                  hp::FECollection<dim>(FE_SimplexP<dim>(d), FE_Q<dim>(d)),
                  cube_and_pyramid);
        test<dim>({1, 0},
                  {1, 0},
                  hp::FECollection<dim>(FE_SimplexP<dim>(d), FE_Q<dim>(d)),
                  cube_and_pyramid);
      }
  }
  deallog.pop();
}
