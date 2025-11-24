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

// Test neighbor_child_on_subface for triangles

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include "../tests.h"


template <int dim>
void
test(const bool refine_triangle, const bool standard_oriented_quad, const bool standard_oriented_tri)
{
  Triangulation<dim> triangulation;
  {
    std::vector<Point<spacedim>> vertices;
    std::vector<CellData<dim>>   cells;
    vertices.push_back(Point<dim>(0, 0));
    vertices.push_back(Point<dim>(1, 0));
    vertices.push_back(Point<dim>(2, 0));
    vertices.push_back(Point<dim>(0, 1));
    vertices.push_back(Point<dim>(1, 1));

    //Quad
    if(standard_oriented_quad)
    {
      CellData<dim> quad;
      quad.vertices = {1,4,3,0};
      cells.push_back(quad);
    }
    else
    {
      CellData<dim> quad;
      quad.vertices = {4,1,0,3};
      cells.push_back(quad);
    }
    //Tri
    if(standard_oriented_tri)
    {
      CellData<dim> tri;
      tri.vertices = {4,1,2};
      cells.push_back(tri);
    }
    else
    {
      CellData<dim> tri;
      tri.vertices = {1,4,2};
      cells.push_back(tri);
    }
    triangulation.create_triangulation(vertices, cells, SubCellData());
  }

  const auto quad = triangulation.begin(0);
  const auto tri = ++ triangulation.begin(0);

  if(refine_triangle)
    tri->set_refine_flag(
      RefinementCase<dim>::isotropic_refinement);
  else
    quad->set_refine_flag(
      RefinementCase<dim>::isotropic_refinement);
  triangulation.execute_coarsening_and_refinement();
 

  const auto cell = refine_triangle ? tri : quad;
  const auto child_1 = cell->face(0)->child(0);
  const auto child_2 = cell->face(0)->child(1);

  const auto neighbor = refine_triangle ? quad : tri;
  auto child_cell_1 = neighbor->neighbor_child_on_subface(0, 0);
  auto child_cell_2 = neighbor->neighbor_child_on_subface(0, 1);

  deallog << child_cell_1->vertex(1) - child_1->vertex(0) << " "
          << child_cell_1->vertex(2) - child_1->vertex(1) << std::endl;

  deallog << child_cell_2->vertex(1) - child_2->vertex(0) << " "
          << child_cell_2->vertex(2) - child_2->vertex(1) << std::endl;
}

int
main()
{
  using namespace dealii;
  initlog();
  test<2>(true, true, true);

  return 0;
}
