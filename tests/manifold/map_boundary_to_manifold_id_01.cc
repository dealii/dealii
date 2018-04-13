//----------------------------  manifold_id_02.cc  ---------------------------
//    Copyright (C) 2018 by the deal.II authors
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  manifold_id_02.cc  ---------------------------


// Test the function map_boundary_to_manifold_ids

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

  // Test setting all manifold ids to 0
  deallog << "All manifold ids to 0" << std::endl;
  auto bids = tria.get_boundary_ids();
  std::vector<types::manifold_id> mids(bids.size(), 0);
  GridTools::map_boundary_to_manifold_ids(bids, mids, tria);
  print_info(tria);

  // Test resetting all boundary ids to 1, and all manifold ids to 2
  mids = std::vector<types::manifold_id>(bids.size(), 2);
  std::vector<types::boundary_id> rbids(bids.size(), 1);
  GridTools::map_boundary_to_manifold_ids(bids, mids, tria, rbids);
  print_info(tria);
}

int main ()
{
  initlog();
  test<1,1>();
  test<1,2>();
  test<1,3>();
  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}

