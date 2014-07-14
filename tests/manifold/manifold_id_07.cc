//----------------------------  manifold_id_05.cc  ---------------------------
//    Copyright (C) 2011, 2013, 2014 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  manifold_id_05.cc  ---------------------------


// Set a manifold id on one of the boundary face, and attach the
// boundary description to it. Refine globally twice and output the mesh.
// Do this building first the Tria, then the manifold, and setting the
// manifold at the end as well, to check proper destruction of the classes.

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

  Point<spacedim> center;
  for (unsigned int i=0; i<dim; ++i)
    center[i] = .5;

  Triangulation<dim,spacedim> tria;
  HyperBallBoundary<dim,spacedim> boundary(center,.5*std::sqrt((double)dim));
  GridGenerator::hyper_cube (tria);
  typename Triangulation<dim,spacedim>::active_cell_iterator cell;

  tria.begin_active()->face(0)->set_manifold_id(1);
  tria.set_manifold(1,boundary);

  tria.refine_global(2);

  for (cell=tria.begin_active(); cell!=tria.end(); ++cell)
    {
      deallog << "C: " << cell << std::endl;
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        deallog << "F: " << cell->face(f) << ", mid: "
                << (int) cell->face(f)->manifold_id() << std::endl;
    }


  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
  tria.set_manifold(1);
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

