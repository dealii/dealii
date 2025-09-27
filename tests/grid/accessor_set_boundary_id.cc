// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// TriaAccessor<0,dim,spacedim>::set_boundary_id() and
// TriaAccessor<0,dim,spacedim>::set_all_boundary_ids() were marked as
// 'const', in contrast to the corresponding functions of faces of
// higher-dimensional objects. This prevented us from writing loops
// over all faces where the 'face' object was marked as 'const'.

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

  std::map<typename Triangulation<dim, spacedim>::face_iterator, int> mymap;

  for (auto &cell : tria.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      {
        face->set_all_boundary_ids(2);
        if (face->at_boundary())
          face->set_boundary_id(1);
      }
  deallog << "OK" << std::endl;
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
