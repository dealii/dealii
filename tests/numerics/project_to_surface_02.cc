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



// test GridTools::project_to_object for quads



#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



class Rotate2d
{
public:
  Rotate2d(const double angle)
    : angle(angle)
  {}

  template <int spacedim>
  Point<spacedim>
  operator()(const Point<spacedim> p) const
  {
    Point<spacedim> q;
    q[0] = std::cos(angle) * p[0] - std::sin(angle) * p[1];
    q[1] = std::sin(angle) * p[0] + std::cos(angle) * p[1];
    for (unsigned d = 2; d < spacedim; ++d)
      q[d] = p[d];
    return q;
  }

private:
  const double angle;
};


template <int dim>
void
do_rotate(Triangulation<dim> &tria)
{
  GridTools::transform(Rotate2d(numbers::PI / 4), tria);
}


void
do_rotate(Triangulation<1> &)
{}



template <int dim>
void
create_triangulation(const bool rotate, Triangulation<dim> &tria)
{
  GridGenerator::hyper_cube(tria, 1., 3.);

  if (rotate)
    do_rotate(tria);
}

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

  for (unsigned int case_no = 0; case_no < 2; ++case_no)
    {
      deallog << "  Case " << case_no << std::endl;
      create_triangulation((case_no == 1), tria);

      const typename Triangulation<dim>::active_cell_iterator cell =
        tria.begin_active();
      Point<dim> trial_point;
      for (unsigned int i = 0; i < dim; ++i)
        trial_point[i] = 1.5;

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
      tria.clear();
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);
  deallog << std::fixed;

  test<2>();
  test<3>();
}
