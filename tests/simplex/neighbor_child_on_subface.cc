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
test(const bool refine_triangle, const bool standard_oriented_quad, const bool standard_oriented_tri)
{
  Triangulation<dim> triangulation;
  {
    std::vector<Point<dim>> vertices;
    std::vector<CellData<dim>>   cells;
    vertices.push_back(Point<dim>(0, 0));
    vertices.push_back(Point<dim>(1, 0));
    vertices.push_back(Point<dim>(2, 0));
    vertices.push_back(Point<dim>(0, 1));
    vertices.push_back(Point<dim>(1, 1));

    if(standard_oriented_quad && standard_oriented_tri)
    {
      CellData<dim> quad;
      quad.vertices = {0,1,4,3};
      cells.push_back(quad);

      CellData<dim> tri;
      tri.vertices = {1,2,4};
      cells.push_back(tri);
    }
    else
    {
      DEAL_II_ASSERT_UNREACHABLE();
    }
    triangulation.create_triangulation(vertices, cells, SubCellData());
  }

  unsigned int quad_face;
  unsigned int tri_face;

  const auto quad = triangulation.begin(0);
  deallog << "quad: " << quad->vertex(0) << ", " << quad->vertex(1) << ", " << quad->vertex(2) << ", " << quad->vertex(3) << std::endl;
  for(auto f : quad->face_indices())
  {
  if(!quad->face(f)->at_boundary())
  {quad_face = f;
  deallog << "quad face " << f << " is not at the boundary" << std::endl; }
   if(quad->face_orientation(f))
   deallog << "quad face " << f << " has standard orientation " << std::endl;
   else
   deallog << "quad face " << f << " has non-standard orientation" << std::endl;
  }

  const auto tri = ++ triangulation.begin(0);
  deallog << "tri: " << tri->vertex(0) << ", " << tri->vertex(1) << ", " << tri->vertex(2) << std::endl;
  for(auto f : tri->face_indices())
  {
    if(!tri->face(f)->at_boundary()){
      tri_face = f;
    
    deallog << "tri face " << f << " is not at the boundary" << std::endl; }
   if(tri->face_orientation(f))
   deallog << "tri face " << f << " has standard orientation " << std::endl;
   else
   deallog << "tri face " << f << " has non-standard orientation" << std::endl;
  }

  if(refine_triangle)
    tri->set_refine_flag(
      RefinementCase<dim>::isotropic_refinement);
  else
    quad->set_refine_flag(
      RefinementCase<dim>::isotropic_refinement);
  triangulation.execute_coarsening_and_refinement();
 

  const auto refined_cell = refine_triangle ? tri : quad;
  const unsigned int face_index = refine_triangle ? tri_face : quad_face;
  const auto child_1 = refined_cell->face(face_index)->child(0);
  const auto child_2 = refined_cell->face(face_index)->child(1);

  const auto neighbor = refine_triangle ? quad : tri;
  auto child_cell_1 = neighbor->neighbor_child_on_subface(face_index, 0);
  auto child_cell_2 = neighbor->neighbor_child_on_subface(face_index, 1);

  deallog << child_cell_1->vertex(0) << " " // << child_1->vertex(0) << " "
          << child_cell_1->vertex(1)  << std::endl; //<< " " << child_1->vertex(1) << std::endl;

  deallog << child_cell_2->vertex(0) << " " // << child_2->vertex(0) << " "
          << child_cell_2->vertex(1) << std::endl; // << " " << child_2->vertex(1) << std::endl;
}

int
main()
{
  using namespace dealii;
  initlog();

  //for(unsigned int i = 0; i < 2; ++i)
  //  for(unsigned int j = 0; j < 2; ++j)
  //    for(unsigned int k = 0; k < 2; ++k)
  //      test<2>(k,i,j);
  test<2>(true,true,true);


  return 0;
}
