// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check VectorTools::project for Vector<double> arguments


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



// define a function that is zero within the square
//   [-sqrt(2)/2,sqrt(2)/2]^d
// which is the domain produced by GridGenerator::hyper_sphere
// when using a linear mapping. we will want to see a solution
// that is not zero when using a higher order mapping and a Q2
// element
template <int dim>
class F : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int = 0) const
  {
    // compute the linfty norm of p
    double m = 0;
    for (unsigned int d = 0; d < dim; ++d)
      m = std::max(m, std::fabs(p[d]));

    // let the value increase linearly away from the square
    return std::max(0., m - std::sqrt(2.) / 2);
  }
};


template <int dim>
void
test()
{
  SphericalManifold<dim> boundary;

  // create 2 triangulations with the
  // same coarse grid, and refine
  // them differently
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  tria.set_manifold(0, boundary);


  FE_Q<dim>       fe(1);
  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe);

  Vector<double> v(dh.n_dofs());

  AffineConstraints<double> cm;
  cm.close();

  // use the implicit Q1 mapping. this will yield a zero solution
  {
    VectorTools::project(dh, cm, QGauss<dim>(3), F<dim>(), v);
    deallog << v.l2_norm() << std::endl;
    Assert(v.l2_norm() == 0, ExcInternalError());
  }

  // use an explicit Q1 mapping. this will yield a zero solution
  {
    VectorTools::project(MappingQ<dim>(1), dh, cm, QGauss<dim>(3), F<dim>(), v);
    deallog << v.l2_norm() << std::endl;
    Assert(v.l2_norm() == 0, ExcInternalError());
  }

  // use an explicit Q2 mapping. this will yield a nonzero solution if it's a
  // straight projection since some of the quadrature points are lying outside
  // the area where the function is zero
  {
    VectorTools::project(MappingQ<dim>(2), dh, cm, QGauss<dim>(3), F<dim>(), v);
    deallog << v.l2_norm() << std::endl;
    Assert(v.l2_norm() != 0, ExcInternalError());
  }

  // use an explicit Q2 mapping but enforce zero boundary values. this will
  // yield a nonzero solution since some of the quadrature points are lying
  // outside the area where the function is zero, even though the values of
  // the DoFs at the boundary are in fact zero (they are interpolated only at
  // points where the function is zero)
  {
    VectorTools::project(
      MappingQ<dim>(2), dh, cm, QGauss<dim>(3), F<dim>(), v, true);
    deallog << v.l2_norm() << std::endl;
    Assert(v.l2_norm() != 0, ExcInternalError());
  }

  // use an explicit Q2 mapping and project to the boundary first. this will
  // yield a nonzero solution since some of the quadrature points are lying
  // outside the area where the function is zero. furthermore, the values on
  // the boundary will be nonzero since we project, even though the
  // *interpolation* of boundary values onto the trace of the Q1 space is zero
  {
    VectorTools::project(MappingQ<dim>(2),
                         dh,
                         cm,
                         QGauss<dim>(3),
                         F<dim>(),
                         v,
                         false,
                         QGauss<dim - 1>(2),
                         true);
    deallog << v.l2_norm() << std::endl;
    Assert(v.l2_norm() != 0, ExcInternalError());
    for (typename DoFHandler<dim>::active_cell_iterator cell =
           dh.begin_active();
         cell != dh.end();
         ++cell)
      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        deallog << cell->vertex(i) << ' ' << v(cell->vertex_dof_index(i, 0))
                << std::endl;
  }


  // same as above, but use a projection with a QTrapezoid formula. this happens
  // to evaluate the function only at points where it is zero, and
  // consequently the values at the boundary should be zero
  {
    VectorTools::project(MappingQ<dim>(2),
                         dh,
                         cm,
                         QGauss<dim>(3),
                         F<dim>(),
                         v,
                         false,
                         QTrapezoid<dim - 1>(),
                         true);
    deallog << v.l2_norm() << std::endl;
    Assert(v.l2_norm() != 0, ExcInternalError());
    for (typename DoFHandler<dim>::active_cell_iterator cell =
           dh.begin_active();
         cell != dh.end();
         ++cell)
      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        deallog << cell->vertex(i) << ' ' << v(cell->vertex_dof_index(i, 0))
                << std::endl;
  }
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
}
