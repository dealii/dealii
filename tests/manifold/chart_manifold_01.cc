// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test that the flat manifold does what it should, where the
// flat manifold is implemented as a ChartManifold with identity
// pull-back and push-forward

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>


template <int dim, int spacedim>
class MyFlatManifold : public ChartManifold<dim,spacedim,spacedim>
{
public:
  virtual std::unique_ptr<Manifold<dim,spacedim> >
  clone () const override
  {
    return std::unique_ptr<Manifold<dim,spacedim> >(new MyFlatManifold());
  }

  virtual
  Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const override
  {
    return space_point;
  }


  virtual
  Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const override
  {
    return chart_point;
  }

  virtual
  DerivativeForm<1,spacedim,spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const override
  {
    DerivativeForm<1,spacedim,spacedim> x;
    for (unsigned int d=0; d<spacedim; ++d)
      x[d][d] = 1;
    return x;
  }
};

// Helper function
template <int dim, int spacedim>
void
test(unsigned int ref=1)
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;

  MyFlatManifold<dim,spacedim> flat_manifold;
  Triangulation<dim,spacedim> tria;
  tria.set_manifold (0, flat_manifold);
  GridGenerator::hyper_cube (tria);
  tria.refine_global(1);

  typename Triangulation<dim,spacedim>::active_cell_iterator
  cell;

  for (cell=tria.begin_active(); cell!=tria.end(); ++cell)
    {
      cell->set_all_manifold_ids (0);

      // check that FlatManifold returns the middle of the cell.
      deallog << "Cell: " << cell << std::endl;
      if (cell->get_manifold().get_new_point_on_cell(cell).distance(cell->center()) > 1e-6)
        {
          deallog << "Default manifold: " << cell->get_manifold().get_new_point_on_cell(cell) << std::endl;
          deallog << "Center of cell  : " << cell->center() << std::endl;
        }
      else
        {
          deallog << "OK!" << std::endl;
        }
    }
}

int
main ()
{
  initlog();

  test<1,1>();
  test<1,2>();
  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}
