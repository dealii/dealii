// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this function tests the correctness of the implementation of matrix free
// operations in getting the function values, the function gradients, and the
// function Laplacians on a cartesian mesh (hyper cube). This tests whether
// cartesian meshes are treated correctly. The test case is without any
// constraints

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include "../tests.h"

#include "interpolate_functions_common.h"


template <int dim>
class CompareFunction : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    double value = 1.2 * p[0];
    for (unsigned int d = 1; d < dim; ++d)
      value -= 2.7 * d * p[d];
    return value;
  }
  virtual Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int) const
  {
    Tensor<1, dim> grad;
    grad[0] = 1.2;
    for (unsigned int d = 1; d < dim; ++d)
      grad[d] = -2.7 * d;
    return grad;
  }
  virtual SymmetricTensor<2, dim>
  hessian(const Point<dim> &p, const unsigned int) const
  {
    return SymmetricTensor<2, dim>();
  }
};


template <int dim, int fe_degree>
void
test()
{
  if (fe_degree == 0)
    return;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  {
    FE_Q<dim>       fe(fe_degree);
    DoFHandler<dim> dof(tria);
    dof.distribute_dofs(fe);

    AffineConstraints<double> constraints;
    constraints.close();
    do_test<dim, fe_degree, double>(dof, constraints);
  }
  {
    FE_DGQ<dim>     fe(fe_degree);
    DoFHandler<dim> dof(tria);
    dof.distribute_dofs(fe);

    AffineConstraints<double> constraints;
    constraints.close();
    do_test<dim, fe_degree, double>(dof, constraints);
  }
  deallog << "Test without templates on FEEvaluation" << std::endl;
  {
    FE_DGQ<dim>     fe(fe_degree);
    DoFHandler<dim> dof(tria);
    dof.distribute_dofs(fe);

    AffineConstraints<double> constraints;
    constraints.close();
    do_test<dim, -1, double>(dof, constraints);
  }
}
