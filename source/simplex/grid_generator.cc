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

#include <deal.II/base/config.h>

#include <deal.II/simplex/grid_generator.h>

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_SIMPLEX_SUPPORT

namespace GridGenerator
{
  template <int dim, int spacedim>
  void
  subdivided_hyper_rectangle_with_simplices(
    Triangulation<dim, spacedim> &   tria,
    const std::vector<unsigned int> &repetitions,
    const Point<dim> &               p1,
    const Point<dim> &               p2,
    const bool                       colorize)
  {
    AssertDimension(dim, spacedim);

    AssertThrow(colorize DEAL_II_EQUALS false, ExcNotImplemented());

    std::vector<Point<spacedim>> vertices;
    std::vector<CellData<dim>>   cells;

    if (dim DEAL_II_EQUALS 2)
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
    else if (dim DEAL_II_EQUALS 3)
      {
        // determine cell sizes
        const Point<dim> dx((p2[0] - p1[0]) / repetitions[0],
                            (p2[1] - p1[1]) / repetitions[1],
                            (p2[2] - p1[2]) / repetitions[2]);

        // create vertices
        for (unsigned int k = 0; k <= repetitions[2]; ++k)
          for (unsigned int j = 0; j <= repetitions[1]; ++j)
            for (unsigned int i = 0; i <= repetitions[0]; ++i)
              vertices.push_back(Point<spacedim>(p1[0] + dx[0] * i,
                                                 p1[1] + dx[1] * j,
                                                 p1[2] + dx[2] * k));

        // create cells
        for (unsigned int k = 0; k < repetitions[2]; ++k)
          for (unsigned int j = 0; j < repetitions[1]; ++j)
            for (unsigned int i = 0; i < repetitions[0]; ++i)
              {
                // create reference HEX cell
                std::array<unsigned int, 8> quad{
                  {(k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                     (j + 0) * (repetitions[0] + 1) + i + 0,
                   (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                     (j + 0) * (repetitions[0] + 1) + i + 1,
                   (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                     (j + 1) * (repetitions[0] + 1) + i + 0,
                   (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                     (j + 1) * (repetitions[0] + 1) + i + 1,
                   (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                     (j + 0) * (repetitions[0] + 1) + i + 0,
                   (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                     (j + 0) * (repetitions[0] + 1) + i + 1,
                   (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                     (j + 1) * (repetitions[0] + 1) + i + 0,
                   (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                     (j + 1) * (repetitions[0] + 1) + i + 1}};

                // TET cell 0
                {
                  CellData<dim> cell;
                  if (((i % 2) + (j % 2) + (k % 2)) % 2 DEAL_II_EQUALS 0)
                    cell.vertices = {{quad[0], quad[1], quad[2], quad[4]}};
                  else
                    cell.vertices = {{quad[0], quad[1], quad[3], quad[5]}};

                  cells.push_back(cell);
                }

                // TET cell 1
                {
                  CellData<dim> cell;
                  if (((i % 2) + (j % 2) + (k % 2)) % 2 DEAL_II_EQUALS 0)
                    cell.vertices = {{quad[2], quad[1], quad[3], quad[7]}};
                  else
                    cell.vertices = {{quad[0], quad[3], quad[2], quad[6]}};
                  cells.push_back(cell);
                }

                // TET cell 2
                {
                  CellData<dim> cell;
                  if (((i % 2) + (j % 2) + (k % 2)) % 2 DEAL_II_EQUALS 0)
                    cell.vertices = {{quad[1], quad[4], quad[5], quad[7]}};
                  else
                    cell.vertices = {{quad[0], quad[4], quad[5], quad[6]}};
                  cells.push_back(cell);
                }

                // TET cell 3
                {
                  CellData<dim> cell;
                  if (((i % 2) + (j % 2) + (k % 2)) % 2 DEAL_II_EQUALS 0)
                    cell.vertices = {{quad[2], quad[4], quad[7], quad[6]}};
                  else
                    cell.vertices = {{quad[3], quad[5], quad[7], quad[6]}};
                  cells.push_back(cell);
                }

                // TET cell 4
                {
                  CellData<dim> cell;
                  if (((i % 2) + (j % 2) + (k % 2)) % 2 DEAL_II_EQUALS 0)
                    cell.vertices = {{quad[1], quad[2], quad[4], quad[7]}};
                  else
                    cell.vertices = {{quad[0], quad[3], quad[6], quad[5]}};
                  cells.push_back(cell);
                }
              }
      }
    else
      {
        AssertThrow(colorize DEAL_II_EQUALS false, ExcNotImplemented());
      }

    // actually create triangulation
    tria.create_triangulation(vertices, cells, SubCellData());
  }

  template <int dim, int spacedim>
  void
  subdivided_hyper_cube_with_simplices(Triangulation<dim, spacedim> &tria,
                                       const unsigned int repetitions,
                                       const double       p1,
                                       const double       p2,
                                       const bool         colorize)
  {
    if (dim DEAL_II_EQUALS 2)
      {
        subdivided_hyper_rectangle_with_simplices(
          tria, {{repetitions, repetitions}}, {p1, p1}, {p2, p2}, colorize);
      }
    else if (dim DEAL_II_EQUALS 3)
      {
        subdivided_hyper_rectangle_with_simplices(
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

#endif

// explicit instantiations
#include "grid_generator.inst"

DEAL_II_NAMESPACE_CLOSE
