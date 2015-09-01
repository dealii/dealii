//------------------------------------------------------------
//    Copyright (C) 2015, the deal.II authors
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//------------------------------------------------------------


// Copy from boundary ids to manifold ids

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

template<int dim, int spacedim>
void print_info(Triangulation<dim,spacedim> &tria)
{
  typename Triangulation<dim,spacedim>::active_cell_iterator cell;

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      deallog << "cell: " << cell
              << ", material_id: "
              << (int)cell->material_id()
              << ", manifold_id: "
              << (int)cell->manifold_id() << std::endl;

      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        deallog << "face: " << cell->face(f)
                << ", boundary_id: "
                << (int)cell->face(f)->boundary_id()
                << ", manifold_id: "
                << (int)cell->face(f)->manifold_id() << std::endl;
    }
}


// Helper function
template <int dim, int spacedim>
void test()
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube (tria, 0., 1.);
  tria.refine_global(1);
  tria.begin_active()->set_material_id(1);

  deallog << "Original mesh ==============================" << std::endl;
  print_info(tria);
  GridTools::copy_material_to_manifold_id(tria);
  deallog << "Copied mesh ================================" << std::endl;
  print_info(tria);
  GridTools::copy_material_to_manifold_id(tria, true);
  deallog << "Copied mesh with boundary  =================" << std::endl;
  print_info(tria);
}

int main ()
{
  initlog(true);

  test<1,1>();
  test<1,2>();
  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}

