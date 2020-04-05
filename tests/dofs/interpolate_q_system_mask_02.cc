// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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



// check that VectorTools::interpolate with a component mask works for
// FE_System(FE_Q(p)) elements correctly on an adaptively refined mesh for
// functions of degree q

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"


template <int dim>
class F : public Function<dim>
{
public:
  F(const unsigned int q, const double adj)
    : Function<dim>(3)
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
  triangulation.refine_global(1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  triangulation.refine_global(1);

  for (unsigned int p = 1; p < 6 - dim; ++p)
    {
      FE_Q<dim>       fe1(p);
      FE_Q<dim>       fe2(p + 1);
      FESystem<dim>   fe(fe1, 2, fe2, 1);
      DoFHandler<dim> dof_handler(triangulation);
      dof_handler.distribute_dofs(fe);

      const double                 adj1 = 0.3;
      ComponentSelectFunction<dim> select_mask1(0, 3);
      ComponentMask                mask1(3, false);
      mask1.set(0, true);

      const double                 adj2 = 1.7;
      ComponentSelectFunction<dim> select_mask2(std::make_pair(1, 3), 3);
      ComponentMask                mask2(3, false);
      mask2.set(1, true);
      mask2.set(2, true);

      AffineConstraints<double> constraints;
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      constraints.close();

      Vector<double> interpolant(dof_handler.n_dofs());
      Vector<float>  error(triangulation.n_active_cells());
      for (unsigned int q = 0; q <= p + 2; ++q)
        {
          // interpolate the function
          VectorTools::interpolate(dof_handler,
                                   F<dim>(q, adj1),
                                   interpolant,
                                   mask1);
          VectorTools::interpolate(dof_handler,
                                   F<dim>(q, adj2),
                                   interpolant,
                                   mask2);
          constraints.distribute(interpolant);


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
