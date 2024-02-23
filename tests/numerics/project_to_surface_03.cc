// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like the _02 test, but project to the faces of a distorted cell



#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
typename Triangulation<dim>::quad_iterator
get_quad_iterator(const typename Triangulation<dim>::active_cell_iterator &cell,
                  const unsigned int quad);

template <>
typename Triangulation<2>::quad_iterator
get_quad_iterator<2>(
  const typename Triangulation<2>::active_cell_iterator &cell,
  const unsigned int /*quad*/)
{
  return *reinterpret_cast<const typename Triangulation<2>::quad_iterator *>(
    &cell);
}

template <>
typename Triangulation<3>::quad_iterator
get_quad_iterator<3>(
  const typename Triangulation<3>::active_cell_iterator &cell,
  const unsigned int                                     quad)
{
  return cell->quad(quad);
}


template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> tria;

  GridGenerator::hyper_cube(tria, 0, 1);

  const typename Triangulation<dim>::active_cell_iterator cell =
    tria.begin_active();

  // distort the cell a bit. all
  // faces but face 0 stay planar;
  // face 0 becomes a saddle
  cell->vertex(0)[0] -= 0.25;
  cell->vertex(6)[0] -= 0.25;
  cell->vertex(2)[0] += 0.25;
  cell->vertex(4)[0] += 0.25;


  for (unsigned int point = 0; point < 9; ++point)
    {
      // choose the 8 vertices of the
      // original unit cell as well
      // as the center point
      const Point<dim> trial_point =
        (point < 8 ? GeometryInfo<dim>::unit_cell_vertex(point) :
                     Point<dim>(.5, .5, .5));

      deallog << "Trial point = " << trial_point << std::endl;

      for (unsigned int e = 0; e < GeometryInfo<dim>::quads_per_cell; ++e)
        {
          const typename Triangulation<dim>::quad_iterator quad =
            get_quad_iterator<dim>(cell, e);

          deallog << "    Quad " << e << ", projected point=";

          const Point<dim> p = GridTools::project_to_object(quad, trial_point);
          deallog << p;
          deallog << "  (quad is from ";
          deallog << quad->vertex(0);
          deallog << " to ";
          deallog << quad->vertex(1);
          deallog << " to ";
          deallog << quad->vertex(2);
          deallog << " to ";
          deallog << quad->vertex(3);
          deallog << ')' << std::endl;

          // now make sure that p is
          // indeed closer to
          // trial_point than any of
          // the vertices of the quad
          for (unsigned int v = 0; v < 4; ++v)
            AssertThrow(p.distance(trial_point) <=
                          quad->vertex(v).distance(trial_point),
                        ExcInternalError());
        }
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);
  deallog << std::fixed;

  test<3>();
}
