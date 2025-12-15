// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test neighbor_child_on_subface for triangles and quads

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include "../tests.h"


template <int dim>
void
test(const bool refine_triangle,
     const bool standard_oriented_quad,
     const bool standard_oriented_tri)
{
  if (refine_triangle)
    deallog << "Refine tri" << std::endl;
  else
    deallog << "Refine quad" << std::endl;
  Triangulation<dim> triangulation;
  {
    std::vector<Point<dim>>    vertices;
    std::vector<CellData<dim>> cells;
    vertices.push_back(Point<dim>(0, 0));
    vertices.push_back(Point<dim>(1, 0));
    vertices.push_back(Point<dim>(2, 0));
    vertices.push_back(Point<dim>(0, 1));
    vertices.push_back(Point<dim>(1, 1));

    if (standard_oriented_quad && !standard_oriented_tri)
      {
        deallog << "Standard quad and non standard tri" << std::endl;
        CellData<dim> quad;
        quad.vertices = {0, 1, 3, 4};
        cells.push_back(quad);

        CellData<dim> tri;
        tri.vertices = {1, 2, 4};
        cells.push_back(tri);
      }
    else if (standard_oriented_quad && standard_oriented_tri)
      {
        deallog << "Standard quad and standard tri" << std::endl;

        CellData<dim> quad;
        quad.vertices = {0, 1, 3, 4};
        cells.push_back(quad);

        CellData<dim> tri;
        tri.vertices = {1, 4, 2};
        cells.push_back(tri);
      }
    else if (!standard_oriented_quad && standard_oriented_tri)
      {
        deallog << "Non standard quad and standard tri" << std::endl;

        CellData<dim> tri;
        tri.vertices = {1, 2, 4};
        cells.push_back(tri);

        CellData<dim> quad;
        quad.vertices = {0, 1, 3, 4};
        cells.push_back(quad);
      }
    else
      {
        return;
      }
    triangulation.create_triangulation(vertices, cells, SubCellData());
  }

  unsigned int quad_face;
  unsigned int tri_face;

  const auto quad =
    standard_oriented_quad ? triangulation.begin(0) : ++triangulation.begin(0);
  for (auto f : quad->face_indices())
    {
      if (!quad->face(f)->at_boundary())
        {
          quad_face = f;
        }
    }

  const auto tri =
    standard_oriented_quad ? ++triangulation.begin(0) : triangulation.begin(0);

  for (auto f : tri->face_indices())
    {
      if (!tri->face(f)->at_boundary())
        {
          tri_face = f;
        }
    }

  if (refine_triangle)
    tri->set_refine_flag(RefinementCase<dim>::isotropic_refinement);
  else
    quad->set_refine_flag(RefinementCase<dim>::isotropic_refinement);
  triangulation.execute_coarsening_and_refinement();


  const auto         non_refined_cell = refine_triangle ? quad : tri;
  const unsigned int face_index_neighbor =
    refine_triangle ? quad_face : tri_face;
  auto child_cell_1 =
    non_refined_cell->neighbor_child_on_subface(face_index_neighbor, 0);
  auto child_cell_2 =
    non_refined_cell->neighbor_child_on_subface(face_index_neighbor, 1);

  if (refine_triangle)
    {
      deallog << child_cell_1->vertex(0) << " " << child_cell_1->vertex(1)
              << " " << child_cell_1->vertex(2) << std::endl;

      deallog << child_cell_2->vertex(0) << " " << child_cell_2->vertex(1)
              << " " << child_cell_2->vertex(2) << std::endl;
    }
  else
    {
      deallog << child_cell_1->vertex(0) << " " << child_cell_1->vertex(1)
              << " " << child_cell_1->vertex(2) << " "
              << child_cell_1->vertex(3) << std::endl;

      deallog << child_cell_2->vertex(0) << " " << child_cell_2->vertex(1)
              << " " << child_cell_2->vertex(2) << " "
              << child_cell_2->vertex(3) << std::endl;
    }
  deallog << std::endl;
}

int
main()
{
  using namespace dealii;
  initlog();

  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      for (unsigned int k = 0; k < 2; ++k)
        test<2>(i, j, k);


  return 0;
}
