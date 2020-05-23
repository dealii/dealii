// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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



// check that the mapping throws an exception for the test case in
// distorted_cells_01.cc

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
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
      case 2:
        deallog << "Twisted cell in " << dim << "d" << std::endl;
        std::swap(vertices[1], vertices[0]);
        break;
      default:
        Assert(false, ExcNotImplemented());
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
    }

  Assert(flag == true, ExcInternalError());

  // now build an FEValues object and compute quadrature points on that cell
  FE_Nothing<dim> dummy;
  QGauss<dim>     quadrature(2);
  FEValues<dim>   fe_values(dummy, quadrature, update_JxW_values);
  // should throw an assertion
  try
    {
      fe_values.reinit(coarse_grid.begin());
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  // only twisted cells for FEValues (pinched cells are OK on Gauss points)
  check<1>(2);
  check<2>(2);
  check<3>(2);
}
