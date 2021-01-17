// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

// Test GridTools::get_coarse_mesh_description()

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



void setup_tria(Triangulation<1> &tria)
{
  GridGenerator::subdivided_hyper_rectangle(
    tria, std::vector<unsigned int>(5u, 1), Point<1>(), Point<1>(1.5), true);
  for (auto &cell : tria.active_cell_iterators())
    cell->set_material_id(
      static_cast<types::material_id>(10.0 * (3.0 + cell->center()[0])));
  for (auto face : tria.active_face_iterators())
    if (0.5 < face->center()[0])
      face->set_all_manifold_ids(42);
  for (auto &cell : tria.active_cell_iterators())
    cell->set_manifold_id(
      static_cast<types::material_id>(10.0 * (2.0 + cell->center()[0])));

  for (const auto bid : tria.get_manifold_ids())
    if (bid != numbers::flat_manifold_id)
      tria.set_manifold(bid, FlatManifold<1>());
}



template <int dim>
void
setup_tria(Triangulation<dim> &tria)
{
  static_assert(dim != 1, "not implemented for dim == 1");
  GridGenerator::hyper_ball(tria);
  // Make the boundary/material id check more robust
  for (auto face : tria.active_face_iterators())
    if (face->at_boundary() && 0.0 < face->center()[0])
      face->set_boundary_id(42);
  for (auto &cell : tria.active_cell_iterators())
    if (cell->center() != Point<dim>())
      {
        const double angle = std::atan2(cell->center()[0], cell->center()[1]);
        cell->set_material_id(
          static_cast<types::material_id>(360 + angle * 180.0 / numbers::PI));
      }
  // Make the manifold id check more robust. This relies on the implementation
  // detail that unassigned manifolds are flat.
  for (auto face : tria.active_face_iterators())
    if (-0.1 < face->center()[0])
      face->set_all_manifold_ids(42);
  for (auto &cell : tria.active_cell_iterators())
    if (cell->center() != Point<dim>())
      {
        const double angle = std::atan2(cell->center()[0], cell->center()[1]);
        cell->set_manifold_id(
          static_cast<types::material_id>(361 - angle * 180.0 / numbers::PI));
      }

  for (const auto bid : tria.get_manifold_ids())
    if (bid != numbers::flat_manifold_id)
      tria.set_manifold(bid, FlatManifold<dim>());
}


template <int dim>
void
test()
{
  {
    // Check that we can copy a coarse Triangulation
    Triangulation<dim> tria;
    setup_tria(tria);

    std::vector<Point<dim>>    vertices;
    std::vector<CellData<dim>> cells;
    SubCellData                subcell_data;
    std::tie(vertices, cells, subcell_data) =
      GridTools::get_coarse_mesh_description(tria);

    Triangulation<dim> tria_2;
    tria_2.create_triangulation(vertices, cells, subcell_data);
    Assert(GridTools::have_same_coarse_mesh(tria, tria_2), ExcInternalError());

    for (const auto bid : tria_2.get_manifold_ids())
      if (bid != numbers::flat_manifold_id)
        tria_2.set_manifold(bid, FlatManifold<dim>());

    GridOut grid_out;
    deallog << "Original Triangulation:" << std::endl;
    grid_out.write_vtk(tria, deallog.get_file_stream());
    deallog << "Triangulation constructed from an unrefined Triangulation:"
            << std::endl;
    grid_out.write_vtk(tria_2, deallog.get_file_stream());
  }

  {
    // Check that we can copy a refined Triangulation
    Triangulation<dim> tria;
    setup_tria(tria);

    tria.refine_global(2);
    unsigned int cell_n = 0;
    for (const auto &cell : tria.active_cell_iterators())
      {
        if (cell_n % 3 == 0)
          cell->set_refine_flag();
        else if (cell_n % 5 == 0)
          cell->set_coarsen_flag();
        ++cell_n;
      }
    tria.execute_coarsening_and_refinement();

    std::vector<Point<dim>>    vertices;
    std::vector<CellData<dim>> cells;
    SubCellData                subcell_data;
    std::tie(vertices, cells, subcell_data) =
      GridTools::get_coarse_mesh_description(tria);

    Triangulation<dim> tria_2;
    tria_2.create_triangulation(vertices, cells, subcell_data);
    Assert(GridTools::have_same_coarse_mesh(tria, tria_2), ExcInternalError());

    GridOut grid_out;
    deallog << "Triangulation constructed from a refined Triangulation:"
            << std::endl;
    grid_out.write_vtk(tria_2, deallog.get_file_stream());
  }
}


int
main()
{
  initlog();

  deallog.push("1d");
  test<2>();
  deallog.pop();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
