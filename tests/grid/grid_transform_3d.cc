// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// a test to check GridTransform in 3d. Produces a simplified version of the
// mesh used in Wolfgang's 2006 NIH proposal


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



int
main()
{
  const unsigned int       dim = 3;
  Point<dim>               origin;
  CylindricalManifold<dim> boundary;
  Triangulation<dim>       tria;
  GridGenerator::cylinder(tria, 1, .7);
  tria.set_manifold(0, boundary);
  tria.refine_global(2);

  // build up a map of vertex indices
  // of boundary vertices to the new
  // boundary points
  std::map<unsigned int, Point<dim>> new_points;

  Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                           endc = tria.end();

  for (cell = tria.begin_active(); cell != endc; ++cell)
    if ((cell->center()[0] - .2) * (cell->center()[0] - .2) +
          (cell->center()[2] - cell->center()[1] / 4) *
            (cell->center()[2] - cell->center()[1] / 4) <
        cell->diameter() * cell->diameter())
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  for (cell = tria.begin_active(); cell != endc; ++cell)
    if ((cell->center()[0] - .2) * (cell->center()[0] - .2) +
          (cell->center()[2] - cell->center()[1] / 4) *
            (cell->center()[2] - cell->center()[1] / 4) <
        cell->diameter() * cell->diameter())
      cell->set_refine_flag();
    else
      cell->set_coarsen_flag();
  tria.execute_coarsening_and_refinement();

  {
    GridOut       grid_out;
    std::ofstream out("output_before.vtk");
    out.precision(5);
    out << std::fixed;
    grid_out.write_vtk(tria, out);
  }


  Triangulation<dim>::face_iterator face;
  for (cell = tria.begin_active(); cell != endc; ++cell)
    for (const unsigned int face_no : GeometryInfo<dim>::face_indices())
      {
        face = cell->face(face_no);
        if (face->at_boundary())
          for (unsigned int vertex_no = 0;
               vertex_no < GeometryInfo<dim>::vertices_per_face;
               ++vertex_no)
            {
              const Point<dim> &old_vertex = face->vertex(vertex_no);

              Point<dim> new_vertex(old_vertex[0],
                                    old_vertex[1],
                                    (old_vertex[0] <= 0 ?
                                       old_vertex[2] / 2 :
                                       old_vertex[2] / (2 + 4 * old_vertex[0] *
                                                              old_vertex[0])));

              new_points[face->vertex_index(vertex_no)] = new_vertex;
            }
      }

  {
    GridOut       grid_out;
    std::ofstream out("output_after.vtk");
    out.precision(5);
    out << std::fixed;
    grid_out.write_vtk(tria, out);
  }

  GridTools::laplace_transform<dim>(new_points, tria, nullptr, true);


  GridOut       grid_out;
  std::ofstream out("output");
  out.precision(5);
  out << std::fixed;
  grid_out.write_gnuplot(tria, out);

  tria.clear();

  return 0;
}
