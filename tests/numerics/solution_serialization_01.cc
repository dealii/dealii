// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// This test is used to make sure that FESeries::Fourier/Legendre
// work with non-primitive elements

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test(Triangulation<dim> &tria, const std::vector<BoundingBox<dim, double>> &bbs)
{
  GridGenerator::subdivided_hyper_cube(tria, 2);

  for (unsigned int i = 0; i < bbs.size(); ++i)
    {
      for (const auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            if (bbs[i].point_inside(cell->center()))
              cell->set_refine_flag();
          }

      tria.execute_coarsening_and_refinement();
    }

  for (unsigned int l = 0; l <= bbs.size(); ++l)
    {
      for (const auto cell : tria.cell_iterators_on_level(l))
        deallog << " " << cell->center() << std::endl;

      deallog << std::endl << std::endl;
    }
}


template <int dim>
void
test()
{
  const unsigned int no_refinement = 2;

  // 1) create triangulation by refining a quadrant over and over;
  // this leads to refining cells on coarser levels

  // ... specify the quadrant via a bounding box
  std::pair<Point<dim, double>, Point<dim, double>> points;
  points.first  = {0.0, 0.5};
  points.second = {0.5, 1.0};

  std::vector<BoundingBox<dim, double>> bbs1(no_refinement,
                                             BoundingBox<dim, double>(points));

  // ... run test
  Triangulation<dim> tria1(
    Triangulation<dim>::MeshSmoothing::limit_level_difference_at_vertices);
  test(tria1, bbs1);

  // 2) create triangulation by refining only cells on the
  // finest level

  // ... specify the cells to refine on each level by a boundig box
  std::vector<BoundingBox<dim, double>> bbs2;

  for (unsigned int l = 0; l < no_refinement; ++l)
    {
      std::vector<Point<dim, double>> points;

      for (const auto cell : tria1.cell_iterators_on_level(l + 1))
        for (const auto v : cell->vertex_indices())
          points.push_back(cell->vertex(v));

      bbs2.emplace_back(points);
    }

  // ... run test
  Triangulation<dim> tria2(
    Triangulation<dim>::MeshSmoothing::limit_level_difference_at_vertices);
  test(tria2, bbs2);
}

int
main()
{
  initlog();

  test<2>();
}
