// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// similar to matrix_vector_faces_15 except FE_DGQArbitraryNodes rather than
// FE_DGQ

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"

#include "create_mesh.h"
#include "matrix_vector_faces_common.h"


template <int dim, int fe_degree>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  create_mesh(tria);

  if (dim == 2)
    tria.refine_global(1);
  {
    typename Triangulation<dim>::active_cell_iterator cell =
      tria.begin_active();
    typename Triangulation<dim>::active_cell_iterator endc    = tria.end();
    unsigned int                                      counter = 0;
    for (; cell != endc; ++cell, ++counter)
      if (cell->is_locally_owned() && counter % 3 == 0)
        cell->set_refine_flag();
    tria.execute_coarsening_and_refinement();
  }

  FE_DGQArbitraryNodes<dim> fe(QGauss<1>(fe_degree + 1));
  DoFHandler<dim>           dof_orig(tria);
  DoFHandler<dim>           dof(tria);
  dof_orig.distribute_dofs(fe);
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

  LinearAlgebra::distributed::Vector<double> in, in_orig, out, out_orig;

  // create MatrixFree with DoFHandler in original numbering
  MatrixFree<dim, double>                          mf_data_orig;
  const QGauss<1>                                  quad(fe_degree + 1);
  typename MatrixFree<dim, double>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
  data.tasks_block_size      = 3;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);
  data.initialize_mapping = true;

  mf_data_orig.reinit(mapping, dof_orig, constraints, quad, data);
  mf_data_orig.initialize_dof_vector(in_orig);
  mf_data_orig.initialize_dof_vector(out_orig);

  // create MatrixFree with renumbered degrees of freedom
  MatrixFree<dim, double> mf_data;
  data.initialize_mapping = false;

  mf_data.reinit(mapping, dof, constraints, quad, data);

  std::vector<types::global_dof_index> renumbering;
  mf_data.renumber_dofs(renumbering);
  dof.renumber_dofs(renumbering);

  data.initialize_mapping = true;
  mf_data.reinit(mapping, dof, constraints, quad, data);

  mf_data.initialize_dof_vector(in);
  mf_data.initialize_dof_vector(out);

  // Set random seed for reproducibility
  Testing::srand(42);
  for (unsigned int i = 0; i < in_orig.locally_owned_size(); ++i)
    {
      const double entry       = Testing::rand() / (double)RAND_MAX;
      in_orig.local_element(i) = entry;
      in(renumbering[i])       = entry;
    }

  MatrixFreeAdvection<dim,
                      fe_degree,
                      fe_degree + 1,
                      double,
                      LinearAlgebra::distributed::Vector<double>>
    mf(mf_data_orig);
  mf.vmult(out_orig, in_orig);

  MatrixFreeAdvection<dim,
                      fe_degree,
                      fe_degree + 1,
                      double,
                      LinearAlgebra::distributed::Vector<double>>
    mf2(mf_data);
  mf2.vmult(out, in);

  for (unsigned int i = 0; i < out.locally_owned_size(); ++i)
    out(renumbering[i]) -= out_orig.local_element(i);

  double diff_norm = out.linfty_norm() / out_orig.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << ' ';

  // test again, now doing matrix-vector product twice
  mf2.vmult(out, in);
  mf2.vmult(out, in);
  for (unsigned int i = 0; i < out.locally_owned_size(); ++i)
    out(renumbering[i]) -= out_orig.local_element(i);
  diff_norm = out.linfty_norm() / out_orig.linfty_norm();
  deallog << diff_norm << std::endl;
}
