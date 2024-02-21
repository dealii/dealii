// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check GeometryInfo::alternating_form_at_vertices for the faces of cells

#include <deal.II/base/geometry_info.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "Checking in " << dim << 'd' << std::endl;

  // check the determinant of the
  // transformation for the reference
  // cell. the determinant should be one in
  // that case
  {
    Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
    for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
      vertices[v] = GeometryInfo<dim>::unit_cell_vertex(v);

    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        Point<dim> face_vertices[GeometryInfo<dim>::vertices_per_face];
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          face_vertices[v] =
            vertices[GeometryInfo<dim>::face_to_cell_vertices(f, v)];

        Tensor<1, dim> alternating_forms[GeometryInfo<dim>::vertices_per_face];
        GeometryInfo<dim - 1>::alternating_form_at_vertices(face_vertices,
                                                            alternating_forms);
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          {
            deallog << "Reference cell: face " << f << ": "
                    << alternating_forms[v] << std::endl;
            AssertThrow(alternating_forms[v].norm() == 1, ExcInternalError());
          }
      }
  }

  // try the same, but move squash the cell
  // in the x-direction by a factor of 10
  {
    Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
    for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
      {
        vertices[v] = GeometryInfo<dim>::unit_cell_vertex(v);
        vertices[v][0] /= 10;
      }

    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        Point<dim> face_vertices[GeometryInfo<dim>::vertices_per_face];
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          face_vertices[v] =
            vertices[GeometryInfo<dim>::face_to_cell_vertices(f, v)];
        Tensor<1, dim> alternating_forms[GeometryInfo<dim>::vertices_per_face];
        GeometryInfo<dim - 1>::alternating_form_at_vertices(face_vertices,
                                                            alternating_forms);
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          {
            deallog << "Squashed cell: face " << f << ": "
                    << alternating_forms[v] << std::endl;
            // faces 0,1 should be
            // unaffected, but all
            // other faces are
            // squashed
            if (f < 2)
              {
                AssertThrow(alternating_forms[v].norm() == 1,
                            ExcInternalError());
              }
            else
              {
                AssertThrow(alternating_forms[v].norm() == 0.1,
                            ExcInternalError());
              }
          }
      }
  }


  // try the same, but move squash the cell
  // in the x-direction by a factor of 10 and
  // rotate it around the z-axis (unless in
  // 1d)
  {
    Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
    for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
      {
        vertices[v] = GeometryInfo<dim>::unit_cell_vertex(v);
        vertices[v][0] /= 10;

        if (dim > 1)
          {
            std::swap(vertices[v][0], vertices[v][1]);
            vertices[v][1] *= -1;
          }
      }

    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        Point<dim> face_vertices[GeometryInfo<dim>::vertices_per_face];
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          face_vertices[v] =
            vertices[GeometryInfo<dim>::face_to_cell_vertices(f, v)];
        Tensor<1, dim> alternating_forms[GeometryInfo<dim>::vertices_per_face];
        GeometryInfo<dim - 1>::alternating_form_at_vertices(face_vertices,
                                                            alternating_forms);
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          {
            deallog << "Squashed+rotated cell: face " << f << ": "
                    << alternating_forms[v] << std::endl;

            // in 2d and 3d, faces 0,1 should be unaffected (just like for the
            // squashed cell, the rotation has nothing to do with face numbers
            // though the direction of the alternating form vector would have
            // rotated along)
            if (f < 2)
              {
                AssertThrow(alternating_forms[v].norm() == 1,
                            ExcInternalError());
              }
            else
              {
                AssertThrow(alternating_forms[v].norm() == 0.1,
                            ExcInternalError());
              }
          }
      }
  }

  // pinched cell
  {
    Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
    for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
      vertices[v] = GeometryInfo<dim>::unit_cell_vertex(v);
    vertices[1] /= 10;

    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        Point<dim> face_vertices[GeometryInfo<dim>::vertices_per_face];
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          face_vertices[v] =
            vertices[GeometryInfo<dim>::face_to_cell_vertices(f, v)];
        Tensor<1, dim> alternating_forms[GeometryInfo<dim>::vertices_per_face];
        GeometryInfo<dim - 1>::alternating_form_at_vertices(face_vertices,
                                                            alternating_forms);
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          deallog << "Pinched cell: face " << f << ": " << alternating_forms[v]
                  << std::endl;
      }
  }


  // inverted cell
  {
    Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
    for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
      vertices[v] = GeometryInfo<dim>::unit_cell_vertex(v);
    std::swap(vertices[0], vertices[1]);

    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        Point<dim> face_vertices[GeometryInfo<dim>::vertices_per_face];
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          face_vertices[v] =
            vertices[GeometryInfo<dim>::face_to_cell_vertices(f, v)];

        Tensor<1, dim> alternating_forms[GeometryInfo<dim>::vertices_per_face];
        GeometryInfo<dim - 1>::alternating_form_at_vertices(face_vertices,
                                                            alternating_forms);
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          deallog << "Inverted cell: face " << f << ": " << alternating_forms[v]
                  << std::endl;
      }
  }
}



int
main()
{
  initlog();

  test<2>();
  test<3>();

  return 0;
}
