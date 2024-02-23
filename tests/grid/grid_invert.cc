// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <string>

#include "../tests.h"

#include "../test_grids.h"

template <int dim>
void
test(bool second_case = false)
{
  std::vector<Point<dim>> vertices(GeometryInfo<dim>::vertices_per_cell);
  vertices[1][1] = 1;
  vertices[2][0] = 1;
  vertices[2][1] = 1;
  vertices[3][0] = 1;
  if (dim == 3)
    {
      for (unsigned int i = 4; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        vertices[i][2] = -1;
      vertices[5][1] = 1;
      vertices[6][0] = 1;
      vertices[6][1] = 1;
      vertices[7][0] = 1;
    }
  std::vector<CellData<dim>> cells(1);
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    cells[0].vertices[i] = i;

  if (dim == 3 && second_case)
    {
      std::swap(cells[0].vertices[1], cells[0].vertices[3]);
      std::swap(cells[0].vertices[5], cells[0].vertices[7]);
      for (unsigned int i = 4; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        vertices[i][2] = 1;
    }

  SubCellData subcelldata;

  TestGrids::reorder_old_to_new_style(cells);
  GridTools::invert_all_negative_measure_cells(vertices, cells);

  Triangulation<dim> tria;
  tria.create_triangulation(vertices, cells, subcelldata);

  std::ostream &logfile = deallog.get_file_stream();
  logfile << "---------------------------------------------" << std::endl
          << "dim=" << dim << (second_case ? ", second case" : ", first case")
          << std::endl
          << std::endl;

  GridOut grid_out;
  grid_out.set_flags(GridOutFlags::Ucd(true));
  grid_out.write_ucd(tria, logfile);
}

int
main()
{
  initlog(false, std::ios_base::fmtflags());
  test<2>();
  test<3>(false);
  test<3>(true);
}
