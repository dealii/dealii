// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
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


// Test that an overlap in boundary and manifold IDs does not lead to
// the boundary being treated as a (curved) manifold upon refinement.

#include "../tests.h"


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <fstream>

using namespace dealii;

template <typename Stream, int dim>
void print_triangulation_data(Stream &stream,
                              const Triangulation<dim> &triangulation)
{
// Boundary id count
  std::map<int,int> boundary_id_count;
  std::map<int,int> manifold_id_count;
  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell!=triangulation.end(); ++cell)
    {
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (cell->face(face)->at_boundary())
            {
              boundary_id_count[cell->face(face)->boundary_id()]++;
              manifold_id_count[cell->face(face)->manifold_id()]++;
            }
          else
            {
              // Prevent double-accounting when face is shared
              if (cell->neighbor_level(face) == cell->level())
                {
                  if (cell->id() < cell->neighbor(face)->id())
                    manifold_id_count[cell->face(face)->manifold_id()]++;
                }
              else
                manifold_id_count[cell->face(face)->manifold_id()]++;

            }
        }
    }

  for (typename std::map<int,int>::const_iterator
       it = boundary_id_count.begin();
       it != boundary_id_count.end(); ++it)
    {
      stream
          << "  Boundary: " << it->first
          << "  ;"
          << "  Number of faces: " << it->second
          << "\n";
    }
  for (typename std::map<int,int>::const_iterator
       it = manifold_id_count.begin();
       it != manifold_id_count.end(); ++it)
    {
      stream
          << "  Manifold: " << it->first
          << "  ;"
          << "  Number of faces: " << it->second
          << "\n";
    }

  stream << std::endl;
}

int main (int argc, char *argv[])
{
  initlog();

  const int dim = 2;
  const unsigned int n_global_refinements = 1;

  // Create a geometry with flat and curved boundaries
  const double inner_radius = 0.25;
  const double outer_radius = 0.5;
  const double tol = 1e-6;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube_with_cylindrical_hole (tria,inner_radius,outer_radius);

  // Enumerate the flat boundaries and the curved one
  // separately. Also provide a manifold ID to the curved
  // surface. Note that the manifold ID coincides with
  // one of the boundary IDs. This is the root cause of the
  // issue...
  const types::manifold_id curved_manifold_id = 1;
  for (typename Triangulation<dim>::active_cell_iterator
       cell = tria.begin_active();
       cell != tria.end(); ++cell)
    {
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (cell->face(face)->at_boundary())
            {
              const Point<dim> pt_centre = cell->face(face)->center();

              for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
                {
                  const Point<dim> pt_vertex = cell->face(face)->vertex(vertex);

                  if (std::abs(pt_vertex.norm() - inner_radius) < tol)
                    {
                      cell->face(face)->set_manifold_id(curved_manifold_id);
                      cell->face(face)->set_boundary_id(5);
                    }
                  else
                    {
                      if (std::abs(pt_centre[0] - outer_radius) < tol)
                        cell->face(face)->set_boundary_id(1);
                      else if (std::abs(pt_centre[0] + outer_radius) < tol)
                        cell->face(face)->set_boundary_id(2);
                      else if (std::abs(pt_centre[1] - outer_radius) < tol)
                        cell->face(face)->set_boundary_id(3);
                      else if (std::abs(pt_centre[1] + outer_radius) < tol)
                        cell->face(face)->set_boundary_id(4);
                    }
                }
            }
        }
    }

  static const SphericalManifold<dim> sphere;
  tria.set_manifold(curved_manifold_id, sphere);
  tria.refine_global(n_global_refinements);

  deallog << "Boundary / manifold information" << std::endl;
  print_triangulation_data(deallog, tria);

  deallog << "Output grid: " << std::endl;
  std::ofstream file_gridout ("grid.vtu");
  GridOut().write_vtu(tria, file_gridout);
  GridOut().write_vtk(tria, deallog.get_file_stream());

  deallog << "OK" << std::endl;

  return 0;
}
