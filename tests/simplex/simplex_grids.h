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


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_description.h>

namespace dealii
{
  namespace GridGenerator
  {
    template <int dim, int spacedim>
    void
    subdivided_hyper_rectangle_with_wedges(
      Triangulation<dim, spacedim>    &tria,
      const std::vector<unsigned int> &repetitions,
      const Point<dim>                &p1,
      const Point<dim>                &p2,
      const bool                       colorize = false)
    {
      AssertDimension(dim, spacedim);

      AssertThrow(colorize == false, ExcNotImplemented());

      std::vector<Point<spacedim>> vertices;
      std::vector<CellData<dim>>   cells;

      if (dim == 3)
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


                  // TRI cell 0
                  {
                    CellData<dim> tri;
                    tri.vertices = {
                      quad[0], quad[1], quad[2], quad[4], quad[5], quad[6]};
                    cells.push_back(tri);
                  }

                  // TRI cell 1
                  {
                    CellData<dim> tri;
                    tri.vertices = {
                      quad[3], quad[2], quad[1], quad[7], quad[6], quad[5]};
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
    subdivided_hyper_cube_with_wedges(Triangulation<dim, spacedim> &tria,
                                      const unsigned int            repetitions,
                                      const double                  p1 = 0.0,
                                      const double                  p2 = 1.0,
                                      const bool colorize              = false)
    {
      if (dim == 3)
        {
          subdivided_hyper_rectangle_with_wedges(
            tria,
            {{repetitions, repetitions, repetitions}},
            {p1, p1, p1},
            {p2, p2, p2},
            colorize);
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }
    }



    template <int dim, int spacedim>
    void
    subdivided_hyper_rectangle_with_pyramids(
      Triangulation<dim, spacedim>    &tria,
      const std::vector<unsigned int> &repetitions,
      const Point<dim>                &p1,
      const Point<dim>                &p2,
      const bool                       colorize = false)
    {
      AssertDimension(dim, spacedim);

      AssertThrow(colorize == false, ExcNotImplemented());

      std::vector<Point<spacedim>> vertices;
      std::vector<CellData<dim>>   cells;

      if (dim == 3)
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
          for (unsigned int k = 0; k < repetitions[2]; ++k)
            for (unsigned int j = 0; j < repetitions[1]; ++j)
              for (unsigned int i = 0; i < repetitions[0]; ++i)
                vertices.push_back(Point<spacedim>(p1[0] + dx[0] * (i + 0.5),
                                                   p1[1] + dx[1] * (j + 0.5),
                                                   p1[2] + dx[2] * (k + 0.5)));

          // create cells
          for (unsigned int k = 0; k < repetitions[2]; ++k)
            for (unsigned int j = 0; j < repetitions[1]; ++j)
              for (unsigned int i = 0; i < repetitions[0]; ++i)
                {
                  // create reference HEX cell
                  std::array<unsigned int, 9> quad{
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
                       (j + 1) * (repetitions[0] + 1) + i + 1,
                     (repetitions[2] + 1) * (repetitions[1] + 1) *
                         (repetitions[0] + 1) +
                       (repetitions[1] * repetitions[0] * k +
                        repetitions[0] * j + i)}};


                  // TRI cell 0
                  {
                    CellData<dim> tri;
                    tri.vertices = {
                      quad[0], quad[2], quad[4], quad[6], quad[8]};
                    cells.push_back(tri);
                  }
                  // TRI cell 1
                  {
                    CellData<dim> tri;
                    tri.vertices = {
                      quad[2], quad[3], quad[6], quad[7], quad[8]};
                    cells.push_back(tri);
                  }
                  // TRI cell 2
                  {
                    CellData<dim> tri;
                    tri.vertices = {
                      quad[6], quad[7], quad[4], quad[5], quad[8]};
                    cells.push_back(tri);
                  }
                  // TRI cell 3
                  {
                    CellData<dim> tri;
                    tri.vertices = {
                      quad[1], quad[5], quad[3], quad[7], quad[8]};
                    cells.push_back(tri);
                  }
                  // TRI cell 4
                  {
                    CellData<dim> tri;
                    tri.vertices = {
                      quad[0], quad[1], quad[2], quad[3], quad[8]};
                    cells.push_back(tri);
                  }
                  // TRI cell 5
                  {
                    CellData<dim> tri;
                    tri.vertices = {
                      quad[0], quad[4], quad[1], quad[5], quad[8]};
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
    subdivided_hyper_cube_with_pyramids(Triangulation<dim, spacedim> &tria,
                                        const unsigned int repetitions,
                                        const double       p1       = 0.0,
                                        const double       p2       = 1.0,
                                        const bool         colorize = false)
    {
      if (dim == 3)
        {
          subdivided_hyper_rectangle_with_pyramids(
            tria,
            {{repetitions, repetitions, repetitions}},
            {p1, p1, p1},
            {p2, p2, p2},
            colorize);
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }
    }



    template <int dim, int spacedim>
    void
    subdivided_hyper_rectangle_with_simplices_mix(
      Triangulation<dim, spacedim>    &tria,
      const std::vector<unsigned int> &repetitions,
      const Point<dim>                &p1,
      const Point<dim>                &p2,
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
          AssertThrow(false, ExcNotImplemented());
        }
    }



    /**
     * Adds a simplex cell to face @p face_no of a hyper_cube cell.
     */
    template <int dim, int spacedim>
    void
    cube_and_pyramid(Triangulation<dim, spacedim> &tria,
                     const unsigned int            face_no = 1)
    {
      Assert(face_no % 2 == 1,
             ExcMessage("Only works for odd face numbers. "
                        "GridReordering::reorder_cells() is not prepared for "
                        "simplices yet (uses GeometryInfo)."));

      // create cube
      Triangulation<dim, spacedim> tria_cube;
      GridGenerator::hyper_cube(tria_cube);
      const auto cube = tria_cube.begin_active();
      AssertIndexRange(face_no, cube->n_faces());

      // extract vertices from specified face, store their midpoint
      std::vector<Point<spacedim>> vertices;
      Point<spacedim>              midpoint;
      const auto                   shared_face = cube->face(face_no);
      for (unsigned int v = 0; v < shared_face->n_vertices(); ++v)
        {
          const auto &vertex = shared_face->vertex(v);
          vertices.push_back(vertex);
          midpoint += vertex;
        }
      midpoint /= vertices.size();

      // add simplex cell in coordinate direction
      const unsigned int coordinate = face_no / 2;
      if constexpr (running_in_debug_mode())
        {
          // all vertices should be in a plane
          for (const auto &vertex : vertices)
            Assert(midpoint[coordinate] == vertex[coordinate],
                   ExcInternalError());
        }

      // add another vertex as tip of triangle/pyramid
      Point<spacedim> tip = midpoint;
      tip[coordinate] += (face_no % 2 == 1) ? 1. : -1.;
      vertices.push_back(tip);

      CellData<dim> simplex(vertices.size());
      std::iota(simplex.vertices.begin(), simplex.vertices.end(), 0);

      Triangulation<dim, spacedim> tria_simplex;
      tria_simplex.create_triangulation(vertices, {simplex}, SubCellData());

      GridGenerator::merge_triangulations(tria_cube, tria_simplex, tria);
    }
  } // namespace GridGenerator
} // namespace dealii
