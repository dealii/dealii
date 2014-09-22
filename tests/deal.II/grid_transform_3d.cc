// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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


// a test to check GridTransform in 3d. Produces a simplified version of the
// mesh used in Wolfgang's 2006 NIH proposal


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>


int main ()
{
  const unsigned int dim=3;
  Point<dim> origin;
  CylinderBoundary<dim> boundary;
  Triangulation<dim> tria;
  tria.set_boundary(0, boundary);
  GridGenerator::cylinder(tria, 1, .7);
  tria.refine_global(2);

  // build up a map of vertex indices
  // of boundary vertices to the new
  // boundary points
  std::map<unsigned int,Point<dim> > new_points;

  Triangulation<dim>::active_cell_iterator cell=tria.begin_active(),
                                           endc=tria.end();

  for (cell=tria.begin_active(); cell!=endc; ++cell)
    if ((cell->center()[0]-.2)*(cell->center()[0]-.2) +
        (cell->center()[2]-cell->center()[1]/4)*(cell->center()[2]-cell->center()[1]/4)
        < cell->diameter() * cell->diameter())
      cell->set_refine_flag ();
  tria.execute_coarsening_and_refinement();

  for (cell=tria.begin_active(); cell!=endc; ++cell)
    if ((cell->center()[0]-.2)*(cell->center()[0]-.2) +
        (cell->center()[2]-cell->center()[1]/4)*(cell->center()[2]-cell->center()[1]/4)
        < cell->diameter() * cell->diameter())
      cell->set_refine_flag ();
    else
      cell->set_coarsen_flag ();
  tria.execute_coarsening_and_refinement();


  Triangulation<dim>::face_iterator face;
  for (cell=tria.begin_active(); cell!=endc; ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        face=cell->face(face_no);
        if (face->at_boundary())
          for (unsigned int vertex_no=0;
               vertex_no<GeometryInfo<dim>::vertices_per_face; ++vertex_no)
            {
              const Point<dim> &old_vertex=face->vertex(vertex_no);

              Point<dim> new_vertex (old_vertex[0],
                                     old_vertex[1],
                                     (old_vertex[0] <= 0 ?
                                      old_vertex[2]/2 :
                                      old_vertex[2]/(2+4*old_vertex[0]*old_vertex[0])));

              new_points[face->vertex_index(vertex_no)] = new_vertex;
            }
      }

  GridGenerator::laplace_transformation (tria, new_points);


  GridOut grid_out;
  std::ofstream out("output");
  out.precision (5);
  out << std::fixed;
  grid_out.write_gnuplot(tria, out);

  tria.clear();

  return 0;
}
