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



// check that we can fix up faces if we get distorted cells because a
// face is out of whack

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
class MyManifold : public Manifold<dim>
{
  virtual std::unique_ptr<Manifold<dim>>
  clone() const override
  {
    return std::unique_ptr<Manifold<dim>>(new MyManifold<dim>());
  }

  virtual Point<dim>
  get_new_point_on_line(
    const typename Triangulation<dim>::line_iterator &line) const override
  {
    deallog << "Finding point between " << line->vertex(0) << " and "
            << line->vertex(1) << std::endl;

    return Point<dim>(0, 0.5, 0.9);
  }

  virtual Point<dim>
  get_new_point_on_quad(
    const typename Triangulation<dim>::quad_iterator &) const override
  {
    Assert(false, ExcInternalError());
    return Point<dim>(0, 0, 1.25);
  }
};



template <int dim>
void
check()
{
  MyManifold<dim> my_manifold;

  // create two cubes
  Triangulation<dim> coarse_grid(Triangulation<dim>::none, true);

  std::vector<unsigned int> sub(dim, 1);
  sub[0] = 2;
  Point<dim> p1(-1, 0, 0), p2(1, 1, 1);
  GridGenerator::subdivided_hyper_rectangle(coarse_grid, sub, p1, p2, true);

  // set bottom middle edge to use MyManifold
  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    for (unsigned int e = 0; e < GeometryInfo<dim - 1>::faces_per_cell; ++e)
      if (coarse_grid.begin_active()->face(f)->line(e)->center()[0] == 0)
        if (coarse_grid.begin_active()->face(f)->line(e)->center()[1] == 0.5)
          if (coarse_grid.begin_active()->face(f)->line(e)->center()[2] == 0)
            coarse_grid.begin_active()->face(f)->line(e)->set_manifold_id(99);
  coarse_grid.set_manifold(99, my_manifold);

  // now try to refine this one
  // cell. we should not get an exception, but keep it to make sure the
  // program still runs
  try
    {
      coarse_grid.refine_global(1);
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
      deallog << "Found " << dcv.distorted_cells.size() << " distorted cells"
              << std::endl;

      Assert(dcv.distorted_cells.size() == 2, ExcInternalError());

      typename Triangulation<dim>::DistortedCellList subset =
        GridTools::fix_up_distorted_child_cells(dcv, coarse_grid);
      deallog << "Found " << subset.distorted_cells.size()
              << " cells that are still distorted" << std::endl;

      Assert(subset.distorted_cells.size() == 0, ExcInternalError());
    }

  Assert(coarse_grid.n_levels() == 2, ExcInternalError());
  Assert(coarse_grid.n_active_cells() == 2 * 1 << dim, ExcInternalError());

  // output the coordinates of the
  // child cells
  GridOut().write_gnuplot(coarse_grid, deallog.get_file_stream());
}


int
main()
{
  initlog();

  check<3>();
}
