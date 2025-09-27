// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that we can fix up faces if we get distorted cells because a
// face is out of whack

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
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
