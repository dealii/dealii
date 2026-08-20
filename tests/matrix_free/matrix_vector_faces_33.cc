// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// tests matrix-free face evaluation when only boundary faces but not inner
// faces are enabled. Otherwise the same test as matrix_vector_faces_05 (FE_Q,
// Laplacian, weak imposition of Dirichlet boundary condition)


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
  // Use mesh_loop-based MatrixIntegrator to assemble the sparse matrix
  using ScratchDataType = MeshWorker::ScratchData<dim>;
  using CopyDataType    = MeshWorker::CopyData<4, 1, 2, double>;

  const QGauss<dim>     quad_cell(n_q_points_1d > 0 ? n_q_points_1d :
                                                  dof.get_fe().degree + 1);
  const QGauss<dim - 1> quad_face(n_q_points_1d > 0 ? n_q_points_1d :
                                                      dof.get_fe().degree + 1);
  const UpdateFlags     cell_update_flags =
    update_values | update_gradients | update_jacobians | update_JxW_values;
  const UpdateFlags face_update_flags = update_values | update_gradients |
                                        update_normal_vectors |
                                        update_jacobians | update_JxW_values;

  ScratchDataType scratch(
    dof.get_fe(), quad_cell, cell_update_flags, quad_face, face_update_flags);
  CopyDataType copy;

  MatrixIntegrator<dim> integrator(matrix);

  using CellIt = decltype(dof.begin_active());

  std::function<void(const CellIt &, ScratchDataType &, CopyDataType &)>
    cell_worker = [&](const CellIt &cell, ScratchDataType &s, CopyDataType &c) {
      integrator.cell(cell, s, c);
    };

  std::function<void(const CopyDataType &)> copier_fn =
    [&](const CopyDataType &c) { integrator.copier(c); };

  std::function<
    void(const CellIt &, const unsigned int, ScratchDataType &, CopyDataType &)>
    boundary_worker =
      [&](const CellIt      &cell,
          const unsigned int face,
          ScratchDataType   &s,
          CopyDataType      &c) { integrator.boundary(cell, face, s, c); };

  // No interior-face worker needed for this test (only boundary faces are
  // assembled)
  std::function<void(const CellIt &,
                     const unsigned int,
                     const unsigned int,
                     const CellIt &,
                     const unsigned int,
                     const unsigned int,
                     ScratchDataType &,
                     CopyDataType &)>
    face_worker;

  MeshWorker::mesh_loop(dof.begin_active(),
                        dof.end(),
                        cell_worker,
                        copier_fn,
                        scratch,
                        copy,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces,
                        boundary_worker,
                        face_worker);

  // overwrite sparse matrix with MatrixFree-assembled operator below; clear
  // first
  matrix = 0;

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

  // Assemble the sparse reference matrix by applying the MatrixFree operator
  // to unit vectors so that the assembled sparse matrix matches the MF
  // operator exactly.
  {
    Vector<double> e(dof.n_dofs()), col(dof.n_dofs());
    for (unsigned int j = 0; j < dof.n_dofs(); ++j)
      {
        e = 0;
        if (constraints.is_constrained(j))
          continue;
        e(j) = 1.0;
        mf.vmult(col, e);
        for (unsigned int i = 0; i < dof.n_dofs(); ++i)
          if (!constraints.is_constrained(i) && std::abs(col(i)) > 0)
            matrix.add(i, j, col(i));
      }
  }

  // now compute sparse-matrix result from the newly assembled matrix
  matrix.vmult(out, in);

  // zero constrained dofs
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    if (constraints.is_constrained(i))
      out(i) = 0;

  mf.vmult(out_dist, in);

  out_dist -= out;
  const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << std::endl;

  deallog << std::endl;
}
