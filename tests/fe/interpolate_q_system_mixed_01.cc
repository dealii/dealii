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



// check that VectorTools::interpolate with a component mask works for
// FE_System(FE_Q(p)) elements correctly on a uniformly refined mesh for
// functions of degree q when a non-interpolating (FE_DGP(p)) element is also
// included in the system.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

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
  F(const unsigned int q, const double adj)
    : Function<dim>(4)
    , q(q)
    , adj(adj)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &v) const
  {
    for (unsigned int c = 0; c < v.size(); ++c)
      {
        v(c) = 0;
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int i = 0; i <= q; ++i)
            v(c) += (d + 1) * (i + 1) * std::pow(p[d], 1. * i) + c + adj;
      }
  }

private:
  const unsigned int q;
  const double       adj;
};



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  for (unsigned int p = 1; p < 6 - dim; ++p)
    {
      FE_Q<dim>       fe_1(p);
      FE_Q<dim>       fe_2(p + 1);
      FE_DGP<dim>     fe_3(p);
      FESystem<dim>   fe(fe_1, 2, fe_2, 1, fe_3, 1);
      DoFHandler<dim> dof_handler(triangulation);
      dof_handler.distribute_dofs(fe);

      // Use constant offset to distinguish between masks
      const double                 adj1 = 0.3;
      ComponentSelectFunction<dim> select_mask1(0, 4);
      ComponentMask                mask1(4, false);
      mask1.set(0, true);

      const double                 adj2 = 1.7;
      ComponentSelectFunction<dim> select_mask2(std::make_pair(1, 3), 4);
      ComponentMask                mask2(4, false);
      mask2.set(1, true);
      mask2.set(2, true);

      ComponentMask mask_fail(4, false);
      mask_fail.set(3, true);

      Vector<double> interpolant(dof_handler.n_dofs());
      Vector<float>  error(triangulation.n_active_cells());
      for (unsigned int q = 0; q <= p + 2; ++q)
        {
          // interpolate the function with mask 1
          VectorTools::interpolate(dof_handler,
                                   F<dim>(q, adj1),
                                   interpolant,
                                   mask1);

          // interpolate the function with mask 2
          VectorTools::interpolate(dof_handler,
                                   F<dim>(q, adj2),
                                   interpolant,
                                   mask2);

          // then compute the interpolation error for mask 1
          VectorTools::integrate_difference(dof_handler,
                                            interpolant,
                                            F<dim>(q, adj1),
                                            error,
                                            QGauss<dim>(q + 2),
                                            VectorTools::L2_norm,
                                            &select_mask1);
          if (q <= p)
            Assert(error.l2_norm() < 1e-12 * interpolant.l2_norm(),
                   ExcInternalError());

          deallog << fe.get_name() << ", Mask 1, P_" << q
                  << ", rel. error=" << error.l2_norm() / interpolant.l2_norm()
                  << std::endl;

          // then compute the interpolation error for mask 2
          VectorTools::integrate_difference(dof_handler,
                                            interpolant,
                                            F<dim>(q, adj2),
                                            error,
                                            QGauss<dim>(q + 2),
                                            VectorTools::L2_norm,
                                            &select_mask2);
          if (q <= p)
            Assert(error.l2_norm() < 1e-12 * interpolant.l2_norm(),
                   ExcInternalError());

          deallog << fe.get_name() << ", Mask 2, P_" << q
                  << ", rel. error=" << error.l2_norm() / interpolant.l2_norm()
                  << std::endl;
        }

      // Test for correct failure
      deallog << fe.get_name() << ", Fails for including non-interpolating FE ";
      try
        {
          VectorTools::interpolate(dof_handler,
                                   F<dim>(0, 0.0),
                                   interpolant,
                                   mask_fail);
        }
      catch (const ExceptionBase &e)
        {
          deallog << "OK" << std::endl;
          deallog << "\tFails with " << e.get_exc_name() << std::endl;
        }
      deallog << std::endl;
    }
}



int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();
  deallog << std::setprecision(3);

  test<1>();
  test<2>();
  test<3>();
}
