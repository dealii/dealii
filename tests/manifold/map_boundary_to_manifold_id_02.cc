// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test the function map_boundary_to_manifold_ids for edges in 3d

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

template<int dim, int spacedim>
void print_info(const Triangulation<dim,spacedim> &tria)
{
  for (auto cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      deallog << "C: " << cell
              << ", manifold id: " << (int)cell->manifold_id() << std::endl;
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        {
          deallog << "f: " << cell->face(f)
                  << ", boundary id: " << (int)cell->face(f)->boundary_id()
                  << ", manifold id: " << (int)cell->face(f)->manifold_id() << std::endl;
          if (dim >=3)
            for (signed int e=0; e<static_cast<signed int>(GeometryInfo<dim>::lines_per_face); ++e)
              deallog << "e: " << cell->face(f)->line(e)
                      << ", boundary id: " << (int)cell->face(f)->line(e)->boundary_id()
                      << ", manifold id: " << (int)cell->face(f)->line(e)->manifold_id() << std::endl;
        }
    }
}

// Helper function
template <int dim, int spacedim>
void test()
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;

  deallog << "Default" << std::endl;
  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube (tria, 0, 1, true);
  print_info(tria);
  deallog << std::endl;

  // Set the edge manifolds to something interesting
  deallog << "Set edge manifold IDs" << std::endl;
  for (auto cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      for (signed int e=0; e<static_cast<signed int>(GeometryInfo<dim>::lines_per_cell); ++e)
        cell->line(e)->set_manifold_id(e);
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        {
          cell->face(f)->set_boundary_id(f);
          for (signed int e=0; e<static_cast<signed int>(GeometryInfo<dim>::lines_per_face); ++e)
            cell->face(f)->line(e)->set_boundary_id(f);
        }
    }
  print_info(tria);
  deallog << std::endl;

  // Test

  // Test setting all manifold ids to an offset of the boundary id
  deallog << "All manifold ids to offset boundary id" << std::endl;
  auto bids = tria.get_boundary_ids();
  std::vector<types::manifold_id> mids(bids.size(), 0);
  for (unsigned int i=0; i<bids.size(); ++i)
    mids[i] = 10 + bids[i];
  GridTools::map_boundary_to_manifold_ids(bids, mids, tria);
  print_info(tria);

  // Test setting all manifold ids to an offset of the boundary id
  // and then resetting the boundary ids
  deallog << "All manifold ids to offset boundary id + boundary id reset" << std::endl;
  for (unsigned int i=0; i<bids.size(); ++i)
    mids[i] = 20 + bids[i];
  std::vector<types::boundary_id> rbids(bids.size(), 1);
  GridTools::map_boundary_to_manifold_ids(bids, mids, tria, rbids);
  print_info(tria);
}

int main ()
{
  initlog();
  test<3,3>();

  return 0;
}

