// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that right number of DoFs is created with hp-constraints for
// FE_SimplexP, FE_WedgeP and FE_PyramidP elements on pure meshes

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test_n_dofs(const hp::FECollection<dim> &fe,
            const Triangulation<dim>    &triangulation,
            const unsigned int           face_no)
{
  DoFHandler<dim> dof_handler(triangulation);
  Assert(fe.size() == 2, ExcInternalError());
  Assert(triangulation.get_reference_cells().size() == 1, ExcInternalError());

  deallog << "Testing number of dofs " << fe[0].get_name() << " vs. "
          << fe[1].get_name() << std::endl;

  const unsigned int n_cells = triangulation.n_active_cells();
  unsigned int       i       = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (i < n_cells / 2)
        {
          cell->set_active_fe_index(0);
        }
      else
        {
          cell->set_active_fe_index(1);
        }
      ++i;
    }

  dof_handler.distribute_dofs(fe);

  unsigned int n_dofs_expected = 0;
  n_dofs_expected += fe[0].get_unit_support_points().size();
  n_dofs_expected += fe[1].get_unit_support_points().size();

  // subtract the interface
  if (fe[0].degree == fe[1].degree)
    {
      // subtract vertices, edges and face
      n_dofs_expected -= fe[0].n_dofs_per_face(face_no);
    }
  else
    {
      // only vertices
      n_dofs_expected -=
        fe[0].reference_cell().face_reference_cell(face_no).n_vertices();
    }

  if (dof_handler.n_dofs() != n_dofs_expected)
    DEAL_II_ASSERT_UNREACHABLE();

  deallog << dof_handler.n_dofs() << " compared to " << n_dofs_expected
          << " dofs" << std::endl;
}

template <int dim>
void
test_mesh(const Triangulation<dim> &tria,
          const unsigned int        face_no,
          const unsigned int        degree1,
          const unsigned int        degree2)
{
  hp::FECollection<dim> fe;

  const auto reference_cells = tria.get_reference_cells();

  if (reference_cells[0].is_simplex())
    {
      fe.push_back(FE_SimplexP<dim>(degree1));
      fe.push_back(FE_SimplexP<dim>(degree2));
    }
  else if (reference_cells[0] == ReferenceCells::Pyramid)
    {
      fe.push_back(FE_PyramidP<dim>(degree1));
      fe.push_back(FE_PyramidP<dim>(degree2));
    }
  else if (reference_cells[0] == ReferenceCells::Wedge)
    {
      fe.push_back(FE_WedgeP<dim>(degree1));
      fe.push_back(FE_WedgeP<dim>(degree2));
    }

  test_n_dofs(fe, tria, face_no);
  deallog << std::endl;
}


int
main()
{
  initlog();

  Triangulation<3>         tria;
  std::vector<Point<3>>    vertices;
  std::vector<CellData<3>> cells;

  {
    // tet mesh with 2 elements
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(1., 1., 0.));
    {
      CellData<3> tet;
      tet.vertices = {0, 1, 2, 3};
      cells.push_back(tet);
    }
    {
      CellData<3> tet;
      tet.vertices = {1, 2, 3, 4};
      cells.push_back(tet);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 4; ++i)
      for (unsigned int j = 1; j < 4; ++j)
        test_mesh<3>(tria, 0, i, j);
  }

  {
    // pyramid mesh with 2 elements
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(0., 0., -1.));
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 1, 2, 3, 4};
      cells.push_back(pyramid);
    }
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 1, 2, 3, 5};
      cells.push_back(pyramid);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 2; ++i)
      for (unsigned int j = 1; j < 2; ++j)
        test_mesh<3>(tria, 0, i, j);
  }

  {
    // pyramid mesh with 2 elements attached at the triangular sides
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(-1., 0., 0.));
    vertices.push_back(Point<3>(-1., -1., 0.));
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 1, 2, 3, 4};
      cells.push_back(pyramid);
    }
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 2, 5, 6, 4};
      cells.push_back(pyramid);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 2; ++i)
      for (unsigned int j = 1; j < 2; ++j)
        test_mesh<3>(tria, 1, i, j);
  }

  {
    // wedge mesh with 2 elements
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(1., 0., 1.));
    vertices.push_back(Point<3>(0., 1., 1.));
    vertices.push_back(Point<3>(1., 1., 1.));
    {
      CellData<3> wedge;
      wedge.vertices = {0, 3, 2, 4, 7, 6};
      cells.push_back(wedge);
    }
    {
      CellData<3> wedge;
      wedge.vertices = {0, 1, 3, 4, 5, 7};
      cells.push_back(wedge);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int j = 1; j < 3; ++j)
      for (unsigned int i = 1; i < 3; ++i)
        test_mesh<3>(tria, 2, i, j);
  }


  {
    // wedge mesh with 2 elements attached at the triangle sides
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(1., 0., 1.));
    vertices.push_back(Point<3>(0., 1., 1.));
    vertices.push_back(Point<3>(0., 0., -1.));
    vertices.push_back(Point<3>(1., 0., -1.));
    vertices.push_back(Point<3>(0., 1., -1.));
    {
      CellData<3> wedge;
      wedge.vertices = {0, 1, 2, 3, 4, 5};
      cells.push_back(wedge);
    }
    {
      CellData<3> wedge;
      wedge.vertices = {0, 1, 2, 6, 7, 8};
      cells.push_back(wedge);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 3; ++i)
      for (unsigned int j = 1; j < 3; ++j)
        test_mesh<3>(tria, 0, i, j);
  }

  return 0;
}
