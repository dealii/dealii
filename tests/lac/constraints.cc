// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2023 by the deal.II authors
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



// create and print a bunch of ConstrainMatrices

#include <deal.II/base/parameter_handler.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"


namespace
{
  /**
   * This file uses a different ordering for the vertices in a hex
   * cell than we usually do in deal.II. The different convention used
   * here originates in what we believed the ordering to be in UCD
   * format, until it was discovered in 2022 that UCD will interpret
   * this ordering to correspond to inverted cells -- as a
   * consequence, the UCD ordering was fixed, but the current file is
   * stuck on the old ordering.
   */
  constexpr std::array<unsigned int, 8> local_vertex_numbering{
    {0, 1, 5, 4, 2, 3, 7, 6}};

  /**
   * Following is a set of functions that reorder the data from the
   * "current" to the "classic" format of vertex numbering of cells
   * and faces. These functions do the reordering of their arguments
   * in-place.
   */
  void
  reorder_old_to_new_style(std::vector<CellData<2>> &cells)
  {
    for (auto &cell : cells)
      std::swap(cell.vertices[2], cell.vertices[3]);
  }


  void
  reorder_old_to_new_style(std::vector<CellData<3>> &cells)
  {
    // undo the ordering above
    unsigned int tmp[GeometryInfo<3>::vertices_per_cell];
    for (auto &cell : cells)
      {
        for (const unsigned int i : GeometryInfo<3>::vertex_indices())
          tmp[i] = cell.vertices[i];
        for (const unsigned int i : GeometryInfo<3>::vertex_indices())
          cell.vertices[local_vertex_numbering[i]] = tmp[i];
      }
  }
} // namespace


void
make_tria(Triangulation<3> &tria, int step)
{
  switch (step)
    {
      case 0:
      case 1:
        {
          // two cells packed behind each
          // other. if step==0, refine back one,
          // otherwise the one in front
          const Point<3> vertices[12]        = {Point<3>(0, 0, 0),
                                         Point<3>(1, 0, 0),
                                         Point<3>(1, 0, 1),
                                         Point<3>(0, 0, 1),

                                         Point<3>(0, 1, 0),
                                         Point<3>(1, 1, 0),
                                         Point<3>(1, 1, 1),
                                         Point<3>(0, 1, 1),

                                         Point<3>(0, 2, 0),
                                         Point<3>(1, 2, 0),
                                         Point<3>(1, 2, 1),
                                         Point<3>(0, 2, 1)};
          const int      cell_vertices[2][8] = {{0, 1, 2, 3, 4, 5, 6, 7},
                                           {4, 5, 6, 7, 8, 9, 10, 11}};
          std::vector<CellData<3>> cells(2, CellData<3>());
          for (unsigned int cell = 0; cell < 2; ++cell)
            for (unsigned int j = 0; j < 8; ++j)
              cells[cell].vertices[j] = cell_vertices[cell][j];
          cells[0].material_id = 0;
          cells[1].material_id = 0;

          reorder_old_to_new_style(cells);
          tria.create_triangulation(std::vector<Point<3>>(&vertices[0],
                                                          &vertices[12]),
                                    cells,
                                    SubCellData()); // no boundary information

          if (step == 0)
            tria.last_active()->set_refine_flag();
          else
            tria.begin_active()->set_refine_flag();
          tria.execute_coarsening_and_refinement();

          break;
        };

      case 2:
      case 3:
        {
          // two cells packed next to each
          // other. if step==2, refine right one,
          // otherwise the left one
          const Point<3> vertices[12]        = {Point<3>(0, 0, 0),
                                         Point<3>(1, 0, 0),
                                         Point<3>(1, 0, 1),
                                         Point<3>(0, 0, 1),

                                         Point<3>(0, 1, 0),
                                         Point<3>(1, 1, 0),
                                         Point<3>(1, 1, 1),
                                         Point<3>(0, 1, 1),

                                         Point<3>(2, 0, 0),
                                         Point<3>(2, 0, 1),
                                         Point<3>(2, 1, 0),
                                         Point<3>(2, 1, 1)};
          const int      cell_vertices[2][8] = {{0, 1, 2, 3, 4, 5, 6, 7},
                                           {1, 8, 9, 2, 5, 10, 11, 6}};
          std::vector<CellData<3>> cells(2, CellData<3>());
          for (unsigned int cell = 0; cell < 2; ++cell)
            for (unsigned int j = 0; j < 8; ++j)
              cells[cell].vertices[j] = cell_vertices[cell][j];
          cells[0].material_id = 0;
          cells[1].material_id = 0;

          reorder_old_to_new_style(cells);
          tria.create_triangulation(std::vector<Point<3>>(&vertices[0],
                                                          &vertices[12]),
                                    cells,
                                    SubCellData()); // no boundary information

          if (step == 2)
            tria.last_active()->set_refine_flag();
          else
            tria.begin_active()->set_refine_flag();
          tria.execute_coarsening_and_refinement();

          break;
        };

      case 4:
      case 5:
        {
          // two cells packed on top of each
          // other. if step==4, refine top one,
          // otherwise the bottom one
          const Point<3> vertices[12]        = {Point<3>(0, 0, 0),
                                         Point<3>(1, 0, 0),
                                         Point<3>(1, 0, 1),
                                         Point<3>(0, 0, 1),

                                         Point<3>(0, 1, 0),
                                         Point<3>(1, 1, 0),
                                         Point<3>(1, 1, 1),
                                         Point<3>(0, 1, 1),

                                         Point<3>(1, 0, 2),
                                         Point<3>(0, 0, 2),
                                         Point<3>(1, 1, 2),
                                         Point<3>(0, 1, 2)};
          const int      cell_vertices[2][8] = {{0, 1, 2, 3, 4, 5, 6, 7},
                                           {3, 2, 8, 9, 7, 6, 10, 11}};
          std::vector<CellData<3>> cells(2, CellData<3>());
          for (unsigned int cell = 0; cell < 2; ++cell)
            for (unsigned int j = 0; j < 8; ++j)
              cells[cell].vertices[j] = cell_vertices[cell][j];
          cells[0].material_id = 0;
          cells[1].material_id = 0;

          reorder_old_to_new_style(cells);
          tria.create_triangulation(std::vector<Point<3>>(&vertices[0],
                                                          &vertices[12]),
                                    cells,
                                    SubCellData()); // no boundary information

          if (step == 4)
            tria.last_active()->set_refine_flag();
          else
            tria.begin_active()->set_refine_flag();
          tria.execute_coarsening_and_refinement();

          break;
        };


      case 6:
      case 7:
      case 8:
        {
          // four cells, with several refined
          // (see below)
          const Point<3> vertices[18] = {Point<3>(0, 0, 0),
                                         Point<3>(1, 0, 0),
                                         Point<3>(1, 0, 1),
                                         Point<3>(0, 0, 1),

                                         Point<3>(0, 1, 0),
                                         Point<3>(1, 1, 0),
                                         Point<3>(1, 1, 1),
                                         Point<3>(0, 1, 1),

                                         Point<3>(2, 0, 0),
                                         Point<3>(2, 0, 1),
                                         Point<3>(2, 1, 0),
                                         Point<3>(2, 1, 1),

                                         Point<3>(0, 2, 0),
                                         Point<3>(1, 2, 0),
                                         Point<3>(1, 2, 1),
                                         Point<3>(0, 2, 1),

                                         Point<3>(2, 2, 0),
                                         Point<3>(2, 2, 1)};

          const int cell_vertices[4][8] = {{0, 1, 2, 3, 4, 5, 6, 7},
                                           {1, 8, 9, 2, 5, 10, 11, 6},
                                           {4, 5, 6, 7, 12, 13, 14, 15},
                                           {5, 10, 11, 6, 13, 16, 17, 14}};
          std::vector<CellData<3>> cells(4, CellData<3>());
          for (unsigned int cell = 0; cell < 4; ++cell)
            for (unsigned int j = 0; j < 8; ++j)
              cells[cell].vertices[j] = cell_vertices[cell][j];
          cells[0].material_id = 0;
          cells[1].material_id = 0;
          cells[2].material_id = 0;
          cells[3].material_id = 0;

          reorder_old_to_new_style(cells);
          tria.create_triangulation(std::vector<Point<3>>(&vertices[0],
                                                          &vertices[18]),
                                    cells,
                                    SubCellData()); // no boundary information

          switch (step)
            {
              case 6:
                tria.begin_active()->set_refine_flag();
                break;

              case 7:
                tria.begin_active()->set_refine_flag();
                (std::next(tria.begin_active()))->set_refine_flag();
                break;
              case 8:
                tria.begin_active()->set_refine_flag();
                (std::next((++(++tria.begin_active()))))->set_refine_flag();
                break;
            };

          tria.execute_coarsening_and_refinement();

          break;
        };
    };
}



int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog.get_file_stream() << std::setprecision(2);

  FiniteElement<3> *fe = nullptr;

  for (unsigned int element = 0; element < 2; ++element)
    {
      switch (element)
        {
          case 0:
            fe = new FE_Q<3>(1);
            break;
          case 1:
            fe = new FE_Q<3>(2);
            break;
        };

      for (int step = 0; step < 9; ++step)
        {
          deallog << "Element=" << element << ", Step=" << step << std::endl;

          Triangulation<3> tria;
          make_tria(tria, step);
          GridOut().write_gnuplot(tria, deallog.get_file_stream());

          DoFHandler<3> dof(tria);
          dof.distribute_dofs(*fe);

          AffineConstraints<double> constraints;
          DoFTools::make_hanging_node_constraints(dof, constraints);
          constraints.close();

          constraints.print(deallog.get_file_stream());

          // release FE
          dof.clear();

          deallog << std::endl;
        };

      delete fe;
    };
}
