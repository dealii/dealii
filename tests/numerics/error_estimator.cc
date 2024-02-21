// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */



#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class MySquareFunction : public Function<dim>
{
public:
  MySquareFunction()
    : Function<dim>(2)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return (component + 1) * p.square();
  }

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const
  {
    values(0) = value(p, 0);
    values(1) = value(p, 1);
  }
};



template <int dim>
Quadrature<dim - 1> &
get_q_face(Function<dim> &)
{
  static QGauss<dim - 1> q(4);
  return q;
}


template <int dim>
void
check()
{
  Functions::CosineFunction<dim> function;

  Triangulation<dim> tr;
  if (dim == 2)
    {
      GridGenerator::hyper_ball(tr, Point<dim>(), 1);
      tr.reset_manifold(0);
    }
  else
    GridGenerator::hyper_cube(tr, -1, 1);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  if (dim == 1)
    tr.refine_global(2);

  FE_Q<dim>       element(QIterated<1>(QTrapezoid<1>(), 3));
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  MappingQ<dim>        mapping(3);
  Quadrature<dim - 1> &q_face = get_q_face(function);

  std::map<types::boundary_id, const Function<dim> *> neumann_bc;
  neumann_bc[0] = &function;

  Vector<double> v(dof.n_dofs());
  VectorTools::interpolate(mapping, dof, function, v);

  Vector<float> error(tr.n_active_cells());

  KellyErrorEstimator<dim>::estimate(
    mapping, dof, q_face, neumann_bc, v, error);

  deallog << "Estimated error:" << std::endl;
  for (unsigned int i = 0; i < error.size(); ++i)
    deallog << error(i) * 100 << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  deallog.push("1d");
  check<1>();
  deallog.pop();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
