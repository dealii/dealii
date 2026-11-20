// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// GridTools::collect_periodic_faces() with a single boundary id relies on
// local face indices of hypercube cells. Check that it reports an error on a
// simplex mesh and that the variant with two boundary ids works there.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();

  Triangulation<2> tria;
  GridGenerator::subdivided_hyper_rectangle_with_simplices(tria,
                                                           {{2, 2}},
                                                           Point<2>(0., 0.),
                                                           Point<2>(1., 1.));

  for (const auto &cell : tria.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        {
          const auto center = face->center();
          if (std::abs(center[0]) < 1e-12)
            face->set_boundary_id(0);
          else if (std::abs(center[0] - 1.) < 1e-12)
            face->set_boundary_id(1);
          else if (std::abs(center[1]) < 1e-12)
            face->set_boundary_id(2);
          else
            face->set_boundary_id(3);
        }

  using FacePairs =
    std::vector<GridTools::PeriodicFacePair<Triangulation<2>::cell_iterator>>;

  for (unsigned int direction = 0; direction < 2; ++direction)
    {
      FacePairs pairs;
      try
        {
          GridTools::collect_periodic_faces(tria, 0, direction, pairs);
          deallog << "direction " << direction << ": no exception" << std::endl;
        }
      catch (const ExceptionBase &)
        {
          deallog << "direction " << direction << ": exception" << std::endl;
        }
    }

  for (unsigned int direction = 0; direction < 2; ++direction)
    {
      FacePairs pairs;
      GridTools::collect_periodic_faces(
        tria, 2 * direction, 2 * direction + 1, direction, pairs);
      deallog << "direction " << direction << ": " << pairs.size() << " pairs"
              << std::endl;
      for (const auto &pair : pairs)
        deallog << pair.cell[0]->face(pair.face_idx[0])->center() << " <-> "
                << pair.cell[1]->face(pair.face_idx[1])->center() << std::endl;
    }
}
