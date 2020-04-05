// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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
