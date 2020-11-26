// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#include <deal.II/grid/grid_generator.h>

namespace dealii
{
  namespace GridGenerator
  {
    template <int dim, int spacedim>
    void
    subdivided_hyper_rectangle_with_simplices_mix(
      Triangulation<dim, spacedim> &   tria,
      const std::vector<unsigned int> &repetitions,
      const Point<dim> &               p1,
      const Point<dim> &               p2,
      const bool                       colorize = false)
    {
      AssertDimension(dim, spacedim);

      AssertThrow(colorize == false, ExcNotImplemented());

      std::vector<Point<spacedim>> vertices;
      std::vector<CellData<dim>>   cells;

      if (dim == 2)
        {
          // determine cell sizes
          const Point<dim> dx((p2[0] - p1[0]) / repetitions[0],
                              (p2[1] - p1[1]) / repetitions[1]);

          // create vertices
          for (unsigned int j = 0; j <= repetitions[1]; ++j)
            for (unsigned int i = 0; i <= repetitions[0]; ++i)
              vertices.push_back(
                Point<spacedim>(p1[0] + dx[0] * i, p1[1] + dx[1] * j));

          // create cells
          for (unsigned int j = 0; j < repetitions[1]; ++j)
            for (unsigned int i = 0; i < repetitions[0]; ++i)
              {
                // create reference QUAD cell
                std::array<unsigned int, 4> quad{{
                  (j + 0) * (repetitions[0] + 1) + i + 0, //
                  (j + 0) * (repetitions[0] + 1) + i + 1, //
                  (j + 1) * (repetitions[0] + 1) + i + 0, //
                  (j + 1) * (repetitions[0] + 1) + i + 1  //
                }};                                       //

                if (j < repetitions[1] / 2 && i < repetitions[0] / 2)
                  {
                    CellData<dim> quad_;
                    quad_.vertices = {quad[0], quad[1], quad[2], quad[3]};
                    cells.push_back(quad_);

                    continue;
                  }

                // TRI cell 0
                {
                  CellData<dim> tri;
                  tri.vertices = {quad[0], quad[1], quad[2]};
                  cells.push_back(tri);
                }

                // TRI cell 1
                {
                  CellData<dim> tri;
                  tri.vertices = {quad[3], quad[2], quad[1]};
                  cells.push_back(tri);
                }
              }
        }
      else
        {
          AssertThrow(colorize == false, ExcNotImplemented());
        }

      // actually create triangulation
      tria.create_triangulation(vertices, cells, SubCellData());
    }


    template <int dim, int spacedim>
    void
    subdivided_hyper_cube_with_simplices_mix(Triangulation<dim, spacedim> &tria,
                                             const unsigned int repetitions,
                                             const double       p1 = 0.0,
                                             const double       p2 = 1.0,
                                             const bool colorize   = false)
    {
      if (dim == 2)
        {
          subdivided_hyper_rectangle_with_simplices_mix(
            tria, {{repetitions, repetitions}}, {p1, p1}, {p2, p2}, colorize);
        }
      else if (dim == 3)
        {
          subdivided_hyper_rectangle_with_simplices_mix(
            tria,
            {{repetitions, repetitions, repetitions}},
            {p1, p1, p1},
            {p2, p2, p2},
            colorize);
        }
      else
        {
          AssertThrow(false, ExcNotImplemented())
        }
    }
  } // namespace GridGenerator
} // namespace dealii
