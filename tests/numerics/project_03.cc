// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


// check VectorTools::project for Vector<double> arguments
// and inhomogeneous constraints


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



// define the multi-linear function x or x*y or x*y*z that we will
// subsequently project onto the ansatz space
template <int dim>
class F : public Function<dim>
{
public:
  F(){};

  virtual double
  value(const Point<dim> &p, const unsigned int = 0) const
  {
    double s = 1;
    for (unsigned int i = 0; i < dim; ++i)
      s *= p[i];
    return s;
  }
};


template <int dim>
void
test()
{
  // create 2 triangulations with the
  // same coarse grid, and refine
  // them differently
  Triangulation<dim> tria;

  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe);

  Vector<double> v(dh.n_dofs());

  const F<dim> approximated_function;

  AffineConstraints<double> cm;
  VectorTools::interpolate_boundary_values(dh, 0, approximated_function, cm);
  cm.close();

  VectorTools::project(dh, cm, QGauss<dim>(3), approximated_function, v);

  for (typename DoFHandler<dim>::active_cell_iterator cell = dh.begin_active();
       cell != dh.end();
       ++cell)
    for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
      {
        deallog << cell->vertex(i) << ' ' << v(cell->vertex_dof_index(i, 0))
                << std::endl;

        // check that the error is
        // somewhat small. it won't
        // be zero since we project
        // and do not interpolate
        if (std::fabs(v(cell->vertex_dof_index(i, 0)) -
                      F<dim>().value(cell->vertex(i))) > 1e-4)
          {
            deallog << "expected value: " << F<dim>().value(cell->vertex(i))
                    << std::endl;
            AssertThrow(false, ExcInternalError());
          }
      }
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
