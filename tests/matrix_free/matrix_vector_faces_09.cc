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



// similar to matrix_vector_faces_03 but we now compare against the combined
// evaluation with gather_evaluate and integrate_scatter

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

#include "create_mesh.h"
#include "matrix_vector_faces_common.h"



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  create_mesh(tria);

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  cell = tria.begin_active();
  endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 0.5)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (fe_degree == 1)
    tria.refine_global(1);
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 9 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      endc                 = tria.end();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (counter % (7 - i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells() << std::endl; std::cout << "Number
  // of degrees of freedom: " << dof.n_dofs() << std::endl; std::cout << "Number
  // of constraints: " << constraints.n_constraints() << std::endl;

  MappingQ<dim> mapping(dof.get_fe().degree + 1);

  Vector<double> in(dof.n_dofs()), out(dof.n_dofs());
  Vector<double> out_dist(out);

  // Set random seed for reproducibility
  Testing::srand(42);
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = Testing::rand() / (double)RAND_MAX;
      in(i)              = entry;
    }

  MatrixFree<dim, double>                          mf_data;
  const QGauss<1>                                  quad(fe_degree + 1);
  typename MatrixFree<dim, double>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
  data.tasks_block_size      = 3;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  mf_data.reinit(mapping, dof, constraints, quad, data);

  MatrixFreeTest<dim, fe_degree, fe_degree + 1> mf(mf_data);
  mf.vmult(out, in);

  MatrixFreeVariant<dim, fe_degree, fe_degree + 1> mf2(mf_data, false);
  mf2.vmult(out_dist, in);
  mf2.vmult(out_dist, in);

  out_dist -= out;
  const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << std::endl;
}
