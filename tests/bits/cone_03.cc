// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2015 by the deal.II authors
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



// check ConeBoundary

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_c1.h>

#include <fstream>



template <int dim>
void check ()
{
  Triangulation<dim> triangulation;
  GridGenerator::truncated_cone (triangulation);
  Point<dim> p1, p2;
  p1[0] = -1;
  p2[0] = 1;
  static const ConeBoundary<dim> boundary (1, 0.5, p1, p2);
  triangulation.set_boundary (0, boundary);

  triangulation.refine_global (2);

  const typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active();
  for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
    {
      typename Triangulation<dim>::face_iterator face=cell->face(face_no);
      typename Boundary<dim>::FaceVertexNormals normals;
      boundary.get_normals_at_vertices(face, normals);
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
        {
          deallog << normals[v] << std::endl;
        }
    }
}


int main ()
{
  initlog();

  check<2> ();
  check<3> ();
}
