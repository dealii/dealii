// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

// check that we can create filtered face iterators from a face iterator
// range. Similar to filtered_iterator_05.

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <set>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  // refine the boundary cells a few times
  for (unsigned int i = 0; i < 5; ++i)
    {
      for (auto &cell : tria.active_cell_iterators())
        cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  // we now have a number of cells, flag them with some manifold ids
  std::size_t i = 0;
  for (auto &cell : tria.active_cell_iterators())
    {
      cell->set_all_manifold_ids(i % 3);
      ++i;
    }

  // Count the faces that are on the boundary and have a manifold_id of 0
  const types::manifold_id manifold_id = 0;
  std::set<typename Triangulation<dim, spacedim>::active_face_iterator>
    face_set;
  for (const auto &cell : tria.active_cell_iterators())
    for (const unsigned int face_n : GeometryInfo<dim>::face_indices())
      if (cell->face(face_n)->at_boundary() &&
          cell->face(face_n)->manifold_id() == manifold_id)
        face_set.insert(cell->face(face_n));

  std::size_t n_filtered_cells = 0;
  for (const auto filtered_cell : filter_iterators(
         tria.active_face_iterators(),
         IteratorFilters::AtBoundary(),
         [=](const typename Triangulation<dim, spacedim>::active_face_iterator
               &face) { return face->manifold_id() == manifold_id; }))
    {
      AssertThrow(face_set.count(filtered_cell) == 1,
                  ExcMessage("Wrong cell filtered."));
      ++n_filtered_cells;
    }
  AssertThrow(n_filtered_cells == face_set.size(),
              ExcMessage("Filtered cells missing."));
}

int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  deallog << "OK" << std::endl;
}
