// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2017 by the deal.II authors
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



// Check periodic boundary faces

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

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
      deallog << periodic_faces[i].cell[0]->index() << " "
              << periodic_faces[i].cell[1]->index() << " "
              << periodic_faces[i].face_idx[0] << " "
              << periodic_faces[i].face_idx[1] << " "
              << periodic_faces[i].orientation << std::endl;
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
