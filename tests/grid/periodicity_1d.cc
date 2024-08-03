// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check periodic boundary faces

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <bitset>
#include <string>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> tr;
  Point<dim>         p0, p1;
  for (unsigned int d = 0; d < dim; ++d)
    p1[d] = 1 + d;
  std::vector<unsigned int> refinements(dim, 3);
  GridGenerator::subdivided_hyper_rectangle(tr, refinements, p0, p1, true);
  tr.refine_global(4 - dim);

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  for (unsigned int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(tr, 2 * d, 2 * d + 1, d, periodic_faces);

  deallog << "Test run in " << dim << " dimensions" << std::endl;
  for (unsigned int i = 0; i < periodic_faces.size(); ++i)
    {
      deallog << periodic_faces[i].cell[0]->index() << ' '
              << periodic_faces[i].cell[1]->index() << ' '
              << periodic_faces[i].face_idx[0] << ' '
              << periodic_faces[i].face_idx[1]
              << ' '
              // maintain old output by converting to a a bitset
              << std::bitset<3>(periodic_faces[i].orientation) << std::endl;
    }
}


int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
}
