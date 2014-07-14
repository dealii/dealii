//----------------------------  manifold_id_02.cc  ---------------------------
//    Copyright (C) 2011, 2013 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  manifold_id_02.cc  ---------------------------


// Test Manifold ID. Now we test the function set_all_manifold_ids(), and verify
// that they are correctly inherited from one cell onward. All manifold ids
// should end up being 1.

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
  for (unsigned int i=0; i<dim; ++i)
    center[i] = .25;

  // const HyperBallBoundary<dim,spacedim> boundary(center,center.norm());
  // triangulation.set_boundary (0, boundary_description);

  typename Triangulation<dim,spacedim>::active_cell_iterator cell;

  tria.begin_active()->set_all_manifold_ids(1);
  tria.refine_global(ref);

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
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1,1>();
  test<1,2>();
  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}

