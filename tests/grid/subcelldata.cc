// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "../test_grids.h"


static unsigned subcells[6][4] = {{0, 1, 2, 3},
                                  {4, 5, 6, 7},
                                  {0, 1, 5, 4},
                                  {1, 5, 6, 2},
                                  {3, 2, 6, 7},
                                  {0, 4, 7, 3}};



template <int dim>
void
test()
{
  Assert(dim == 2 || dim == 3, ExcNotImplemented());

  std::vector<Point<dim>> vertices(GeometryInfo<dim>::vertices_per_cell);
  vertices[0][0] = 0;
  vertices[0][1] = 0;
  vertices[1][0] = 2;
  vertices[1][1] = 1;
  vertices[2][0] = 3;
  vertices[2][1] = 3;
  vertices[3][0] = 0;
  vertices[3][1] = 1;
  if (dim == 3)
    {
      // for the new numbering
      //       for (unsigned int i=0; i<4; ++i)
      //  {
      //    vertices[i+4]=vertices[i];
      //    vertices[i+4](2)=1;
      //  }
      // for the old numbering
      for (unsigned int i = 0; i < 4; ++i)
        {
          std::swap(vertices[i][1], vertices[i][2]);
          vertices[i + 4]    = vertices[i];
          vertices[i + 4][1] = 1;
        }
    }

  std::vector<CellData<dim>> cells(1);
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    cells[0].vertices[i] = i;
  cells[0].material_id = 0;

  SubCellData subcelldata;
  if (dim == 2)
    {
      subcelldata.boundary_lines.resize(GeometryInfo<dim>::faces_per_cell);
      for (const unsigned int i : GeometryInfo<dim>::face_indices())
        {
          subcelldata.boundary_lines[i].vertices[0] = i;
          subcelldata.boundary_lines[i].vertices[1] = (i + 1) % 4;
          subcelldata.boundary_lines[i].material_id = 10 * i + 1;
        }
    }
  else if (dim == 3)
    {
      subcelldata.boundary_quads.resize(GeometryInfo<dim>::faces_per_cell);
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face;
               ++i)
            subcelldata.boundary_quads[f].vertices[i] = subcells[f][i];
          subcelldata.boundary_quads[f].material_id = 10 * f + 1;
        }
    }

  TestGrids::reorder_old_to_new_style(cells);
  Triangulation<dim> tria;
  tria.create_triangulation(vertices, cells, subcelldata);

  GridOutFlags::Ucd ucd_flags(true, true);
  GridOut           grid_out;
  grid_out.set_flags(ucd_flags);
  grid_out.write_ucd(tria, deallog.get_file_stream());

  //   std::ofstream gnuplot_file("subcelldata.gnuplot");
  //   grid_out.write_gnuplot(tria, gnuplot_file);
  //   std::ofstream ucd_file("subcelldata.inp");
  //   grid_out.write_ucd(tria, ucd_file);
}


int
main()
{
  initlog();

  test<2>();
  test<3>();

  return 0;
}
