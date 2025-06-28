// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
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
// function Laplacians on a hyperball mesh for FE_Q_Hierarchical elements
// (no symmetry in transformation). The test case includes hanging node
// constraints, but no constraints from boundary conditions

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "get_functions_common.h"



template <int dim, int fe_degree>
void
test()
{
  using number = double;
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);

  // refine first and last cell
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(4 - dim);

  FE_Q_Hierarchical<dim> fe(fe_degree);
  DoFHandler<dim>        dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells()
  //          << std::endl;
  // std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
  // std::cout << "Number of constraints: " << constraints.n_constraints() <<
  // std::endl;

  Vector<number> solution(dof.n_dofs());

  // create vector with random entries
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = random_value<double>();
      solution(i)        = entry;
    }
  constraints.distribute(solution);

  MatrixFree<dim, number> mf_data;
  deallog << "Test with fe_degree " << fe_degree << std::endl;
  const QGauss<1>                                  quad(fe_degree + 1);
  MappingQ<dim>                                    mapping(4);
  typename MatrixFree<dim, number>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
  data.mapping_update_flags  = update_gradients | update_hessians;
  mf_data.reinit(mapping, dof, constraints, quad, data);
  MatrixFreeTest<dim, fe_degree, fe_degree + 1, number> mf(mf_data, mapping);
  mf.test_functions(solution);
}
