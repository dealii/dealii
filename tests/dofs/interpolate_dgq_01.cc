// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that VectorTools::interpolate works for FE_DGQ(p) elements correctly on
// a uniformly refined mesh for functions of degree q

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"


template <int dim>
class F : public Function<dim>
{
public:
  F(const unsigned int q)
    : q(q)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    double v = 0;
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int i = 0; i <= q; ++i)
        v += (d + 1) * (i + 1) * std::pow(p[d], 1. * i);
    return v;
  }

private:
  const unsigned int q;
};



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  for (unsigned int p = 1; p < 7 - dim; ++p)
    {
      FE_DGQ<dim>     fe(p);
      DoFHandler<dim> dof_handler(triangulation);
      dof_handler.distribute_dofs(fe);

      Vector<double> interpolant(dof_handler.n_dofs());
      Vector<float>  error(triangulation.n_active_cells());
      for (unsigned int q = 0; q <= p + 2; ++q)
        {
          // interpolate the function
          VectorTools::interpolate(dof_handler, F<dim>(q), interpolant);

          // then compute the interpolation error
          VectorTools::integrate_difference(dof_handler,
                                            interpolant,
                                            F<dim>(q),
                                            error,
                                            QGauss<dim>(q + 2),
                                            VectorTools::L2_norm);
          if (q <= p)
            Assert(error.l2_norm() < 1e-12 * interpolant.l2_norm(),
                   ExcInternalError());

          deallog << fe.get_name() << ", P_" << q
                  << ", rel. error=" << error.l2_norm() / interpolant.l2_norm()
                  << std::endl;
        }
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<1>();
  test<2>();
  test<3>();
}
