// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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


// write a file in the VTK format, then read it back in, and check
// that the triangulations are identical.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"

template <int dim, int spacedim>
void
check(Triangulation<dim, spacedim> &tria)
{
  std::ofstream outfile("grid.vtk");
  GridOut       go;
  go.write_vtk(tria, outfile);
  outfile.close();

  Triangulation<dim, spacedim> tria2;
  GridIn<dim, spacedim>        gi;
  std::ifstream                infile("grid.vtk");
  gi.attach_triangulation(tria2);
  gi.read_vtk(infile);

  deallog << "Testing Triangulation<" << dim << "," << spacedim << ">"
          << std::endl;
  AssertDimension(tria.n_vertices(), tria2.n_vertices());
  AssertDimension(tria.n_active_cells(), tria2.n_active_cells());
  auto cell2 = tria2.begin_active();
  for (auto cell1 : tria.active_cell_iterators())
    {
      AssertDimension(cell1->material_id(), cell2->material_id());
      AssertDimension(cell1->manifold_id(), cell2->manifold_id());

      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        {
          auto obj1 = cell1->vertex_index(i);
          auto obj2 = cell2->vertex_index(i);
          AssertDimension(obj1, obj2);

          auto p1 = cell1->vertex(i);
          auto p2 = cell2->vertex(i);
          Assert(p1.distance(p2) == 0, ExcInternalError());
        }
      if (dim == 3)
        for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_cell; ++i)
          {
            auto obj1 = cell1->line(i);
            auto obj2 = cell2->line(i);
            AssertDimension(obj1->manifold_id(), obj2->manifold_id());
            AssertDimension(obj1->boundary_id(), obj2->boundary_id());
          }
      if (dim > 1)
        for (const unsigned int i : GeometryInfo<dim>::face_indices())
          {
            auto obj1 = cell1->face(i);
            auto obj2 = cell2->face(i);
            AssertDimension(obj1->boundary_id(), obj2->boundary_id());
            AssertDimension(obj1->manifold_id(), obj2->manifold_id());
          }
      ++cell2;
    }
  deallog << "OK" << std::endl;
}

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;

  GridGenerator::hyper_cube(tria, 0, 1, true);

  tria.refine_global(1);

  tria.begin_active()->set_all_manifold_ids(5);
  tria.begin_active()->set_material_id(3);

  if (dim > 1)
    tria.begin_active()->face(0)->set_all_manifold_ids(6);

  check(tria);
}

int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
