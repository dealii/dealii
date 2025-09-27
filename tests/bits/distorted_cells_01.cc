// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that indeed Triangulation::create_triangulation throws an
// exception if we have distorted cells

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



// create a (i) pinched cell (where two vertices coincide), or (ii)
// twisted cell (where two vertices are swapped)
template <int dim>
void
check(const unsigned int testcase)
{
  std::vector<Point<dim>> vertices;
  for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
    vertices.push_back(GeometryInfo<dim>::unit_cell_vertex(v));

  switch (testcase)
    {
      case 1:
        deallog << "Pinched cell in " << dim << 'd' << std::endl;
        vertices[0] = vertices[1];
        break;
      case 2:
        deallog << "Twisted cell in " << dim << 'd' << std::endl;
        std::swap(vertices[0], vertices[1]);
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }


  std::vector<CellData<dim>> cells;
  {
    CellData<dim> cell;
    for (const unsigned int j : GeometryInfo<dim>::vertex_indices())
      cell.vertices[j] = j;
    cells.push_back(cell);
  }

  Triangulation<dim> coarse_grid(Triangulation<dim>::none, true);

  bool flag = false;
  try
    {
      coarse_grid.create_triangulation(vertices, cells, SubCellData());
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
      flag = true;

      deallog << dcv.distorted_cells.size() << " distorted cells" << std::endl;
      Assert(dcv.distorted_cells.front() == coarse_grid.begin(0),
             ExcInternalError());
    }
  catch (ExceptionBase &exc)
    {
      deallog << exc.get_exc_name() << std::endl;
      flag = true;
    }

  Assert(flag == true, ExcInternalError());
}


int
main()
{
  initlog();

  for (unsigned int testcase = 1; testcase <= 2; ++testcase)
    {
      check<1>(testcase);
      check<2>(testcase);
      check<3>(testcase);
    }
}
