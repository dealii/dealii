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



// tests matrix-free face evaluation when only boundary faces but not inner
// faces are enabled. Otherwise the same test as matrix_vector_faces_05 (FE_Q,
// Laplacian, weak imposition of Dirichlet boundary condition)

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_q.h>

#include "../tests.h"

#include "matrix_vector_faces_common.h"

template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5 - dim);

  FE_Q<dim>       fe(fe_degree);
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

  constexpr unsigned int n_q_points_1d = fe_degree + 1;

  // assemble sparse matrix with MeshWorker
  SparsityPattern      sparsity;
  SparseMatrix<double> matrix;
  {
    DynamicSparsityPattern d_sparsity(dof.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof, d_sparsity);
    sparsity.copy_from(d_sparsity);
  }
  matrix.reinit(sparsity);
  MeshWorker::IntegrationInfoBox<dim> info_box;
  UpdateFlags                         update_flags =
    update_values | update_gradients | update_jacobians;
  info_box.add_update_flags_all(update_flags);
  info_box.initialize_gauss_quadrature(n_q_points_1d,
                                       n_q_points_1d,
                                       n_q_points_1d);
  info_box.initialize(dof.get_fe(), mapping);

  MeshWorker::DoFInfo<dim> dof_info(dof);

  MeshWorker::Assembler::MatrixSimple<SparseMatrix<double>> assembler;
  assembler.initialize(matrix);

  MatrixIntegrator<dim> integrator;
  MeshWorker::integration_loop<dim, dim>(
    dof.begin_active(), dof.end(), dof_info, info_box, integrator, assembler);

  matrix.vmult(out, in);

  // zero constrained dofs
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    if (constraints.is_constrained(i))
      out(i) = 0;

  MatrixFree<dim, double> mf_data;
  const QGauss<1>         quad(n_q_points_1d > 0 ? n_q_points_1d :
                                           dof.get_fe().degree + 1);
  typename MatrixFree<dim, double>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
  data.tasks_block_size      = 3;
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  mf_data.reinit(mapping, dof, constraints, quad, data);

  // Check that only boundary faces are set up as requested
  Assert(mf_data.n_inner_face_batches() == 0, ExcInternalError());
  Assert(mf_data.n_boundary_face_batches() > 0, ExcInternalError());

  MatrixFreeTest<dim, fe_degree, n_q_points_1d, double, Vector<double>, 1> mf(
    mf_data);
  mf.vmult(out_dist, in);

  out_dist -= out;
  const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << std::endl;

  deallog << std::endl;
}
