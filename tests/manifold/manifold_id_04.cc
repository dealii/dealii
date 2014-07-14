//----------------------------  manifold_id_03.cc  ---------------------------
//    Copyright (C) 2011, 2013 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  manifold_id_03.cc  ---------------------------


// Test Manifold ID. Now we set differently the id of the first cell and its faces, and verify
// that the ids are correctly inherited from one cell onward.

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube (tria);
  Point<spacedim> center;

  typename Triangulation<dim,spacedim>::active_cell_iterator cell;


  tria.begin_active()->set_manifold_id(3);
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    tria.begin_active()->face(f)->set_manifold_id(2);

  tria.refine_global(1);

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      deallog << "C: " << cell
              << ", mid: " << (int)cell->manifold_id() << std::endl;
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        deallog << "f: " << cell->face(f)
                << ", mid: " << (int)cell->face(f)->manifold_id() << std::endl;
    }

}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(1);
  deallog.threshold_double(1.e-10);

  test<1,1>();
  test<1,2>();
  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}

