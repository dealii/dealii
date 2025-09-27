/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2018 - 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Test n_faces, n_active_faces, and n_raw_faces also in 1D.

#include <deal.II/base/logstream.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  tria.begin_active()->set_coarsen_flag();
  (std::next((tria.begin_active())))->set_coarsen_flag();
  tria.execute_coarsening_and_refinement();

  std::map<typename Triangulation<dim, spacedim>::face_iterator, unsigned int>
    mymap;

  unsigned int face_counter = 0;
  for (auto &cell : tria.active_cell_iterators())
    {
      for (const unsigned int i : GeometryInfo<dim>::face_indices())
        {
          auto face = cell->face(i);
          if (mymap.find(cell->face(i)) == mymap.end())
            {
              face_counter++;
              mymap[cell->face(i)] = face->index();
              if (face->index() != static_cast<int>(cell->face_index(i)))
                deallog << "Different indices!" << std::endl;
            }
        }
    }

  deallog << "Counted: " << face_counter
          << ", active: " << tria.n_active_faces()
          << ", total: " << tria.n_faces() << ", raw: " << tria.n_raw_faces()
          << std::endl;
}

int
main()
{
  initlog();
  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
