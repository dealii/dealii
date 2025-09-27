// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// TriaAccessor<0,dim,spacedim> and TriaAccessor<0,1,spacedim> lacked
// get_triangulation() member functions because they are not derived
// from TriaAccessorBase. This led to awkward follow-up errors in
// strange places when using them in certain contexts.

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
  tria.refine_global(1);

  std::map<typename Triangulation<dim, spacedim>::face_iterator, int> mymap;

  for (auto &cell : tria.active_cell_iterators())
    {
      for (auto &face : cell->face_iterators())
        {
          if (mymap.find(face) == mymap.end())
            {
              mymap[face] = cell->index();
            }
        }
    }

  for (auto it : mymap)
    deallog << it.first << ": " << it.second << std::endl;
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
