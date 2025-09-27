// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test GridTools::invert_cells_with_negative_measure() with a mixed mesh.

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <string>

#include "../tests.h"

void
test()
{
  deal_II_exceptions::disable_abort_on_exception();

  // test a grid with a pyramid that cannot be fixed
  const static int dim = 3;
  try
    {
      std::vector<Point<dim>> vertices(10);
      // hypercube
      vertices[0] = Point<dim>(0, 0, 0);
      vertices[1] = Point<dim>(1, 0, 0);
      vertices[2] = Point<dim>(0, 1, 0);
      vertices[3] = Point<dim>(1, 1, 0);
      vertices[4] = Point<dim>(0, 0, 1);
      vertices[5] = Point<dim>(1, 0, 1);
      vertices[6] = Point<dim>(0, 1, 1);
      vertices[7] = Point<dim>(1, 1, 1);
      // pyramid on top
      vertices[8] = Point<dim>(0.75, 0.5, 1.5);
      // pyramid on right - a mess of an element
      vertices[9] = Point<dim>(1.0, 1.25, 0.5);


      std::vector<CellData<dim>> cells(3);
      cells[0].vertices = {0, 1, 2, 3, 4, 5, 6, 7};
      // this one is twisted
      cells[1].vertices = {4, 5, 7, 6, 8};
      // this one is twisted in an unfixable way
      cells[2].vertices = {1, 3, 5, 7, 9};

      SubCellData       subcelldata;
      const std::size_t n_cells_inverted =
        GridTools::invert_cells_with_negative_measure(vertices, cells);

      deallog << "We inverted " << n_cells_inverted << " cell(s)." << std::endl;
    }
  catch (const ExcGridHasInvalidCell &exception)
    {
      deallog << "Successfully failed to set up a mesh with a pyramid "
              << "pointing into another element" << std::endl;
    }

  // test with a fixable pyramid
  {
    std::vector<Point<dim>> vertices(10);
    // hypercube
    vertices[0] = Point<dim>(0, 0, 0);
    vertices[1] = Point<dim>(1, 0, 0);
    vertices[2] = Point<dim>(0, 1, 0);
    vertices[3] = Point<dim>(1, 1, 0);
    vertices[4] = Point<dim>(0, 0, 1);
    vertices[5] = Point<dim>(1, 0, 1);
    vertices[6] = Point<dim>(0, 1, 1);
    vertices[7] = Point<dim>(1, 1, 1);
    // pyramid on top
    vertices[8] = Point<dim>(0.75, 0.5, 1.5);
    // pyramid on right
    vertices[9] = Point<dim>(1.25, 0.75, 0.5);


    std::vector<CellData<dim>> cells(3);
    cells[0].vertices = {0, 1, 2, 3, 4, 5, 6, 7};
    // this one is twisted
    cells[1].vertices = {4, 5, 7, 6, 8};
    // this one is not twisted
    cells[2].vertices = {5, 1, 7, 3, 9};

    SubCellData       subcelldata;
    const std::size_t n_cells_inverted =
      GridTools::invert_cells_with_negative_measure(vertices, cells);

    deallog << "We inverted " << n_cells_inverted << " cell(s)." << std::endl;

    Triangulation<dim> tria;
    tria.create_triangulation(vertices, cells, subcelldata);

    std::ostream &logfile = deallog.get_file_stream();
    logfile << "---------------------------------------------" << std::endl
            << std::endl
            << std::endl;
    GridOut grid_out;
    grid_out.write_vtk(tria, logfile);
  }

  // test with a twisted and untwisted wedge
  {
    std::vector<Point<dim>> vertices(12);
    // hypercube
    vertices[0] = Point<dim>(0, 0, 0);
    vertices[1] = Point<dim>(1, 0, 0);
    vertices[2] = Point<dim>(0, 1, 0);
    vertices[3] = Point<dim>(1, 1, 0);
    vertices[4] = Point<dim>(0, 0, 1);
    vertices[5] = Point<dim>(1, 0, 1);
    vertices[6] = Point<dim>(0, 1, 1);
    vertices[7] = Point<dim>(1, 1, 1);
    // wedge on top
    vertices[8] = Point<dim>(0.0, 0.5, 1.5);
    vertices[9] = Point<dim>(1.0, 0.6, 1.6);
    // wedge on right
    vertices[10] = Point<dim>(1.0, 0.5, 0.0);
    vertices[11] = Point<dim>(1.0, 0.6, 0.9);


    std::vector<CellData<dim>> cells(3);
    cells[0].vertices = {0, 1, 2, 3, 4, 5, 6, 7};
    // this one is twisted
    cells[1].vertices = {5, 7, 9, 4, 6, 8};
    // this one is not twisted
    cells[2].vertices = {1, 10, 3, 5, 11, 7};

    SubCellData       subcelldata;
    const std::size_t n_cells_inverted =
      GridTools::invert_cells_with_negative_measure(vertices, cells);

    deallog << "We inverted " << n_cells_inverted << " cell(s)." << std::endl;

    Triangulation<dim> tria;
    tria.create_triangulation(vertices, cells, subcelldata);

    std::ostream &logfile = deallog.get_file_stream();
    logfile << "---------------------------------------------" << std::endl
            << std::endl
            << std::endl;
    GridOut grid_out;
    grid_out.write_vtk(tria, logfile);
  }
}

int
main()
{
  initlog();
  test();
}
