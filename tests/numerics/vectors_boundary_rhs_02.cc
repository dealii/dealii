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



// check VectorTools::create_boundary_right_hand_side, this time for
// non-primitive elements


#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_raviart_thomas.h>
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
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return (component + 1) * p.square();
  }

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const
  {
    for (unsigned int d = 0; d < dim; ++d)
      values(d) = value(p, d);
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
  FE_RaviartThomas<dim> element(0);
  DoFHandler<dim>       dof(tr);
  dof.distribute_dofs(element);

  // use a more complicated mapping
  // of the domain and a quadrature
  // formula suited to the elements
  // we have here
  MappingQ<dim> mapping(3);

  QGauss<dim - 1> quadrature(3);

  Vector<double> rhs(dof.n_dofs());
  VectorTools::create_boundary_right_hand_side(dof,
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
  deallog << std::setprecision(4);
  deallog << std::fixed;

  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
