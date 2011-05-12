//----------------------------  transform_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  transform_01.cc  ---------------------------

// check GridTools::transform with the mesh used in the "Possibilities for
// extensions" section of step-38. The test exists because the function
// originally did not allow application to meshes that were already refined.


#include "../tests.h"

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>


template <int dim>
Point<dim> warp (const Point<dim> &p)
{
  Point<dim> q = p;
  q[dim-1] *= 10;

  if (dim >= 2)
    q[0] += 2*std::sin(q[dim-1]);
  if (dim >= 3)
    q[1] += 2*std::cos(q[dim-1]);

  return q;
}


template <int dim, int spacedim>
void save_mesh(const Triangulation<dim,spacedim>& tria)
{
  GridOut grid_out;
  grid_out.write_gnuplot (tria, deallog.get_file_stream());
}


int main ()
{
  std::ofstream logfile("transform_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Triangulation<2,3> triangulation;
  
  HyperBallBoundary<3> boundary_description;
  Triangulation<3> volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);
  
  volume_mesh.set_boundary (1, boundary_description);
  volume_mesh.set_boundary (0, boundary_description);
  volume_mesh.refine_global (3);
  
  static HyperBallBoundary<3-1,3> surface_description;
  triangulation.set_boundary (1, surface_description);
  triangulation.set_boundary (0, surface_description);
  
  std::set<unsigned char> boundary_ids;
  boundary_ids.insert(0);
  
  GridTools::extract_boundary_mesh (volume_mesh, triangulation,
				    boundary_ids);
  triangulation.set_boundary (1);
  triangulation.set_boundary (0);
  GridTools::transform (&warp<3>, triangulation);

  deallog << "Surface mesh has " << triangulation.n_active_cells()
	  << " cells."
	  << std::endl;
  save_mesh (triangulation);
  
  return 0;
}
