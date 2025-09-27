// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check VectorTools::create_right_hand_side


#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

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
void
check()
{
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

  // create a system element composed
  // of one Q1 and one Q2 element
  FESystem<dim>   element(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  // use a more complicated mapping
  // of the domain and a quadrature
  // formula suited to the elements
  // we have here
  MappingQ<dim> mapping(3);

  QGauss<dim> quadrature(3);

  Vector<double> rhs(dof.n_dofs());
  VectorTools::create_right_hand_side(dof,
                                      quadrature,
                                      MySquareFunction<dim>(),
                                      rhs);
  for (unsigned int i = 0; i < rhs.size(); ++i)
    deallog << rhs(i) << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(8);
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
