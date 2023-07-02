// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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



#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <string>

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

template <int dim>
void
test(bool second_case = false)
{
  std::vector<Point<dim>> vertices(GeometryInfo<dim>::vertices_per_cell);
  vertices[1](1) = 1;
  vertices[2](0) = 1;
  vertices[2](1) = 1;
  vertices[3](0) = 1;
  if (dim == 3)
    {
      for (unsigned int i = 4; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        vertices[i](2) = -1;
      vertices[5](1) = 1;
      vertices[6](0) = 1;
      vertices[6](1) = 1;
      vertices[7](0) = 1;
    }
  std::vector<CellData<dim>> cells(1);
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    cells[0].vertices[i] = i;

  if (dim == 3 && second_case)
    {
      std::swap(cells[0].vertices[1], cells[0].vertices[3]);
      std::swap(cells[0].vertices[5], cells[0].vertices[7]);
      for (unsigned int i = 4; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        vertices[i](2) = 1;
    }

  SubCellData subcelldata;

  reorder_old_to_new_style(cells);
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
