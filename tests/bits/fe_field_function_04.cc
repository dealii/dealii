// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// FEFieldFunction ran into an assertion after
// Mapping::transform_real_to_unit_cell started throwing exceptions
// when it couldn't find the point on the reference cell that belongs
// to a given point, rather than just silently giving up

#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class F : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    return p.square();
  }
};


template <int dim>
void
test()
{
  const SphericalManifold<dim> boundary_description;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_ball(triangulation);
  triangulation.set_manifold(0, boundary_description);
  triangulation.refine_global(1);

  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, F<dim>(), solution);

  Functions::FEFieldFunction<2> fe_function(dof_handler, solution);
  std::vector<Point<dim>>       points;

  // add a bunch of points. all
  // points are inside the circle.
  // the problem happens because we
  // walk over a bunch of cells in
  // the process of finding all of
  // these points and then realize
  // when we get to the one at the
  // end that the coordinates for
  // this point can't be found in the
  // cell we have touched last (it's
  // too far away from that cell, and
  // the inverse mapping does not
  // converge
  for (unsigned int i = 0; i < 20; ++i)
    for (unsigned int j = 0; j < 20; ++j)
      points.push_back(Point<dim>(-0.7 + i * 0.07, -0.7 + j * 0.07));
  points.push_back(Point<dim>(-0.27999999999999992, -0.62999999999999989));

  std::vector<double> m(points.size());
  fe_function.value_list(points, m);

  for (unsigned int i = 0; i < m.size(); ++i)
    AssertThrow(std::fabs(m[i] - points[i].square()) <
                    1e-10 * std::fabs(m[i] + points[i].square()) ||
                  std::fabs(m[i] + points[i].square()) < 1e-30,
                ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<2>();

  return 0;
}
