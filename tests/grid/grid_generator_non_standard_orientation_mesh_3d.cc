// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/*
 * Small test to analyse the equivalence of the normal component
 * on the element edges for the Raviart-Thomas elements.
 */



#include "../tests.h"

#define PRECISION 8

#include <deal.II/grid/grid_generator.h>

// STL
#include <map>


template <int dim>
void
plot_all_info(const Triangulation<dim> &tria)
{
  for (const auto &cell : tria.active_cell_iterators())
    {
      CellId current_cell_id(cell->id());

      deallog << "CellId = " << current_cell_id << std::endl
              << "   {index -> face_orientation | face_flip | face_rotation}: "
              << std::endl;
      for (unsigned int face_index = 0;
           face_index < GeometryInfo<dim>::faces_per_cell;
           ++face_index)
        {
          deallog << "      {" << face_index << " -> "
                  << cell->face_orientation(face_index) << " | "
                  << cell->face_flip(face_index) << " | "
                  << cell->face_rotation(face_index) << " } " << std::endl;
        } // face_index

      deallog << "   line orientation: {  ";
      for (unsigned int line_index = 0;
           line_index < GeometryInfo<dim>::lines_per_cell;
           ++line_index)
        {
          deallog << (cell->line_orientation(line_index) ==
                      numbers::default_geometric_orientation)
                  << "  ";
        } // line_index
      deallog << '}' << std::endl << std::endl;
    } // cell
}


int
main(int /*argc*/, char ** /*argv*/)
{
  initlog();
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;

  /*
   * 3D test only
   */
  const int          dim = 3;
  Triangulation<dim> tria_test;

  deallog << "Testing 3D mesh for orientation tests:" << std::endl;

  for (unsigned int config_switch = 0; config_switch < 8; ++config_switch)
    {
      tria_test.clear();

      bool face_orientation = (((config_switch / 4) % 2) == 1);
      bool face_flip        = (((config_switch / 2) % 2) == 1);
      bool face_rotation    = ((config_switch % 2) == 1);

      bool manipulate_first_cube = true;

      GridGenerator::non_standard_orientation_mesh(tria_test,
                                                   face_orientation,
                                                   face_flip,
                                                   face_rotation,
                                                   manipulate_first_cube);

      plot_all_info(tria_test);

      deallog << "****************" << std::endl;
    } // config_switch++

  deallog << "Testing 3D orientation mesh done." << std::endl;

  return 0;
}
