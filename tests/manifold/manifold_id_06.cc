//----------------------------  manifold_id_06.cc  ---------------------------
//    Copyright (C) 2011, 2013, 2014 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  manifold_id_06.cc  ---------------------------


// Set a manifold id on the boundary faces of a small cell, and change also
// the interior boundaries.

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>

// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;

  Point<spacedim> center;
  for (unsigned int i=0; i<spacedim; ++i)
    center[i] = .25;

  double radius=center.norm();

  HyperBallBoundary<dim,spacedim> boundary(center, .25*std::sqrt((double)spacedim));
  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube (tria);
  typename Triangulation<dim,spacedim>::active_cell_iterator cell;

  tria.refine_global(1);

  for (cell=tria.begin_active(); cell!=tria.end(); ++cell)
    if (dim<spacedim && cell->center().distance(center)<radius)
      cell->set_all_manifold_ids(1);
    else
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->center().distance(center)< radius)
          cell->face(f)->set_all_manifold_ids(1);

  tria.set_manifold(1,boundary);
  tria.refine_global(2);

  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
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

