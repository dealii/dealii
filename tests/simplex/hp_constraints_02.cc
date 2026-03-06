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


// like hp_constraints_01 but on mixed meshes
// check that right number of DoFs is created with hp-constraints for
// FE_SimplexP, FE_WedgeP and FE_PyramidP elements on mixed meshes

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
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
  Assert(triangulation.get_reference_cells().size() == 2, ExcInternalError());

  deallog << "Testing number of dofs " << fe[0].get_name() << " vs. "
          << fe[1].get_name() << std::endl;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->reference_cell() == fe[0].reference_cell())
        cell->set_active_fe_index(0);
      else if (cell->reference_cell() == fe[1].reference_cell())
        cell->set_active_fe_index(1);
      else
        DEAL_II_ASSERT_UNREACHABLE();
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
test_mixed_mesh(const Triangulation<dim> &tria,
                const unsigned int        face_no,
                const unsigned int        degree,
                const unsigned int        degree_neighbor,
                const bool                swap_fe_indices)
{
  hp::FECollection<dim> fe;
  const auto            reference_cells = tria.get_reference_cells();

  for (unsigned int i = 0; i < reference_cells.size(); ++i)
    {
      const unsigned int current_degree = (i == 0) ? degree : degree_neighbor;
      const unsigned int index =
        swap_fe_indices ? reference_cells.size() - 1 - i : i;
      const auto reference_cell = reference_cells[index];

      if (reference_cell.is_hyper_cube())
        fe.push_back(FE_Q<dim>(current_degree));
      else if (reference_cell.is_simplex())
        fe.push_back(FE_SimplexP<dim>(current_degree));
      else if (reference_cell == ReferenceCells::Pyramid)
        fe.push_back(FE_PyramidP<dim>(current_degree));
      else if (reference_cell == ReferenceCells::Wedge)
        fe.push_back(FE_WedgeP<dim>(current_degree));
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
    // tet + pyramid mesh
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(-1., 0., 0.));

    CellData<3> pyramid;
    pyramid.vertices = {0, 1, 2, 3, 4};
    cells.push_back(pyramid);


    CellData<3> tet;
    tet.vertices = {0, 1, 4, 5};
    cells.push_back(tet);

    tria.create_triangulation(vertices, cells, SubCellData());
    for (unsigned int i = 1; i < 4; ++i)
      {
        test_mixed_mesh<3>(tria, 1, i, 1, false);
        test_mixed_mesh<3>(tria, 1, 1, i, true);
      }
  }

  {
    // tet + wedge mesh
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
    {
      CellData<3> wedge;
      wedge.vertices = {0, 1, 2, 3, 4, 5};
      cells.push_back(wedge);
    }
    {
      CellData<3> tet;
      tet.vertices = {0, 1, 2, 6};
      cells.push_back(tet);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 4; ++i)
      for (unsigned int j = 1; j < 3; ++j)
        {
          test_mixed_mesh<3>(tria, 1, i, j, false);
          test_mixed_mesh<3>(tria, 1, j, i, true);
        }
  }

  {
    // hex + pyramid mesh
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(0., 0., -1.));
    vertices.push_back(Point<3>(1., 0., -1.));
    vertices.push_back(Point<3>(0., 1., -1.));
    vertices.push_back(Point<3>(1., 1., -1.));
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 1, 2, 3, 4};
      cells.push_back(pyramid);
    }
    {
      CellData<3> hex;
      hex.vertices = {0, 1, 2, 3, 5, 6, 7, 8};
      cells.push_back(hex);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 4; ++i)
      {
        test_mixed_mesh<3>(tria, 0, 1, i, false);
        test_mixed_mesh<3>(tria, 0, i, 1, true);
      }
  }

  {
    // hex + wedge mesh
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0.5, 0.5));
    vertices.push_back(Point<3>(1., 0.5, 0.5));
    vertices.push_back(Point<3>(0., 0., -1.));
    vertices.push_back(Point<3>(1., 0., -1.));
    vertices.push_back(Point<3>(0., 1., -1.));
    vertices.push_back(Point<3>(1., 1., -1.));
    {
      CellData<3> wedge;
      wedge.vertices = {0, 4, 2, 1, 5, 3};
      cells.push_back(wedge);
    }
    {
      CellData<3> hex;
      hex.vertices = {0, 1, 2, 3, 6, 7, 8, 9};
      cells.push_back(hex);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 4; ++i)
      for (unsigned int j = 1; j < 3; ++j)
        {
          test_mixed_mesh<3>(tria, 2, j, i, false);
          test_mixed_mesh<3>(tria, 2, i, j, true);
        }
  }

  {
    // pyramid + wedge mesh at quad face
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0.5, 0.5));
    vertices.push_back(Point<3>(1., 0.5, 0.5));
    vertices.push_back(Point<3>(0., 0., -1.));
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 1, 2, 3, 6};
      cells.push_back(pyramid);
    }
    {
      CellData<3> wedge;
      wedge.vertices = {0, 4, 2, 1, 5, 3};
      cells.push_back(wedge);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 3; ++i)
      {
        test_mixed_mesh<3>(tria, 0, 1, i, false);
        test_mixed_mesh<3>(tria, 2, i, 1, true);
      }
  }

  {
    // pyramid + wedge mesh at triangular face
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(-1., 0., 0.));
    vertices.push_back(Point<3>(-1., 1., 0.));
    vertices.push_back(Point<3>(-1., 0., 1.));
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 1, 2, 3, 4};
      cells.push_back(pyramid);
    }
    {
      CellData<3> wedge;
      wedge.vertices = {0, 4, 2, 5, 6, 7};
      cells.push_back(wedge);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 3; ++i)
      {
        test_mixed_mesh<3>(tria, 1, 1, i, false);
        test_mixed_mesh<3>(tria, 1, i, 1, true);
      }
  }

  return 0;
}
