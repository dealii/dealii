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



// The MatrixCreator::create_laplace_matrix function overload that also
// assembles a right hand side vector had a bug in that the statement that
// assembled the rhs vector was nested in the wrong loop. this was fixed by
// Moritz' commit 14428
//
// the function internally has four branches, with different code used
// for the cases with/without coefficient and scalar/vector-valued
// finite element. we test these four cases through the _01, _02, _03,
// and _04 tests. the version without creating a right hand side
// vector is tested in the _0[1234]a tests



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

#include <deal.II/numerics/matrix_tools.h>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> tr;
  if (dim == 2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1, 1);
  tr.reset_manifold(0);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  if (dim == 1)
    tr.refine_global(2);

  hp::FECollection<dim> element;
  element.push_back(FE_Q<dim>(1));
  element.push_back(FE_Q<dim>(2));
  element.push_back(FE_Q<dim>(QIterated<1>(QTrapezoid<1>(), 3)));

  DoFHandler<dim> dof(tr);

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end();
       ++cell)
    cell->set_active_fe_index(Testing::rand() % element.size());

  dof.distribute_dofs(element);

  // use a more complicated mapping
  // of the domain and a quadrature
  // formula suited to the elements
  // we have here
  MappingQ<dim> mapping(3);
  QGauss<dim>   quadrature(6);

  // create sparsity pattern. note
  // that different components should
  // not couple, so use pattern
  SparsityPattern sparsity(dof.n_dofs(), dof.n_dofs());
  DoFTools::make_sparsity_pattern(dof, sparsity);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();
  constraints.condense(sparsity);
  sparsity.compress();

  SparseMatrix<double> matrix;
  matrix.reinit(sparsity);

  Functions::ExpFunction<dim> rhs_function;

  Vector<double> rhs(dof.n_dofs());

  MatrixTools::create_laplace_matrix(hp::MappingCollection<dim>(mapping),
                                     dof,
                                     hp::QCollection<dim>(quadrature),
                                     matrix,
                                     rhs_function,
                                     rhs,
                                     &rhs_function);

  // since we only generate
  // output with two digits after
  // the dot, and since matrix
  // entries are usually in the
  // range of 1 or below,
  // multiply matrix by 100 to
  // make test more sensitive
  deallog << "Matrix: " << std::endl;
  for (SparseMatrix<double>::const_iterator p = matrix.begin();
       p != matrix.end();
       ++p)
    deallog << p->value() * 100 << std::endl;

  deallog << "RHS vector: " << std::endl;
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    deallog << rhs(i) << std::endl;
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
