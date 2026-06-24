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


// like hp_constraints_02 but for FE_SimplexP and FE_Q connected by an edge
// check that right number of DoFs is created with hp-constraints

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test_n_dofs(const hp::FECollection<dim> &fe,
            const Triangulation<dim>    &triangulation,
            const bool                   equidistant_support_points)
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

  // 2 vertices
  n_dofs_expected -= 2;
  // subtract the line
  if (fe[0].degree == fe[1].degree)
    {
      // for degrees larger than 3 check if the support points are equidistant
      if (fe[0].degree > 2 && !equidistant_support_points)
        {
          // nothing to do
        }
      else
        // subtract the edge dofs
        n_dofs_expected -= fe[0].n_dofs_per_line();
    }

  if (dof_handler.n_dofs() != n_dofs_expected)
    DEAL_II_ASSERT_UNREACHABLE();

  deallog << dof_handler.n_dofs() << " compared to " << n_dofs_expected
          << " dofs" << std::endl;
}


template <int dim>
void
test_mixed_mesh(const Triangulation<dim> &tria,
                const unsigned int        degree,
                const unsigned int        degree_neighbor,
                const bool                swap_fe_indices,
                const bool                equidistant_support_points)
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
        {
          if (equidistant_support_points)
            fe.push_back(
              FE_Q<dim>(QIterated<1>(QTrapezoid<1>(), current_degree)));
          else
            fe.push_back(FE_Q<dim>(current_degree));
        }
      else
        fe.push_back(FE_SimplexP<dim>(current_degree));
    }

  test_n_dofs(fe, tria, equidistant_support_points);
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
    vertices.push_back(Point<3>(1., 0., 1.));
    vertices.push_back(Point<3>(0., 1., 1.));
    vertices.push_back(Point<3>(1., 1., 1.));
    vertices.push_back(Point<3>(-0.5, 0.5, 0.5));
    vertices.push_back(Point<3>(-0.5, 0.5, 1.5));


    CellData<3> hex;
    hex.vertices = {0, 1, 2, 3, 4, 5, 6, 7};
    cells.push_back(hex);


    CellData<3> tet;
    tet.vertices = {4, 6, 8, 9};
    cells.push_back(tet);

    tria.create_triangulation(vertices, cells, SubCellData());
    for (unsigned int i = 1; i < 4; ++i)
      for (unsigned int j = 1; j < 4; ++j)
        {
          test_mixed_mesh<3>(tria, i, j, false, false);
          test_mixed_mesh<3>(tria, i, j, false, true);
          test_mixed_mesh<3>(tria, j, i, true, false);
          test_mixed_mesh<3>(tria, j, i, true, true);
        }
  }

  return 0;
}
