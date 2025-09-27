// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test for flag copy_boundary_ids in merge_triangulations()

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
print_boundary_ids(const dealii::Triangulation<dim> &tria)
{
  for (const auto &face : tria.active_face_iterators())
    if (face->at_boundary())
      deallog << "Face at " << face->center()
              << " has id: " << face->boundary_id() << "\n";
  deallog << "-----------------" << std::endl;
}


int
main()
{
  initlog();

  const unsigned int dim = 2;

  Triangulation<dim> tria1;
  Triangulation<dim> tria2;
  Triangulation<dim> tria3;
  Triangulation<dim> tria;
  GridGenerator::hyper_rectangle(tria1, {-1., -1.}, {0., 1.}, true);
  GridGenerator::hyper_rectangle(tria2, {0., -1.}, {1., 1.}, true);
  for (auto &face : tria2.active_face_iterators())
    if (face->at_boundary())
      face->set_boundary_id(face->boundary_id() + dim * 2);

  deallog << "Initial IDs of Triangulation 1" << std::endl;
  print_boundary_ids(tria1);
  deallog << "Initial IDs of Triangulation 2" << std::endl;
  print_boundary_ids(tria2);

  deallog << "Merge two rectangles, one face is eliminated" << std::endl;
  GridGenerator::merge_triangulations(tria1, tria2, tria, 1e-6, false, true);
  print_boundary_ids(tria);

  deallog
    << "Merge two rectangles, due to hanging node there is no internal face"
    << std::endl;
  tria.clear();
  tria2.refine_global(1);
  GridGenerator::flatten_triangulation(tria2, tria3);
  GridGenerator::merge_triangulations(tria1, tria3, tria, 1e-6, false, true);
  print_boundary_ids(tria);

  deallog << std::endl;
}
