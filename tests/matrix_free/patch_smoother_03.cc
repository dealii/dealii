// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// The same as patch_smoother_02, but using Dynamic PatchDistributor.
// check FEPatchEvaluation by comparing A*x computed locally with A*x computed
// globally via FEEvaluation

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_patch_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/patch_distributors.h>
#include <deal.II/matrix_free/patch_storage.h>

#include "../tests.h"


template <int dim, int fe_degree>
void
test()
{
  deallog << "Running test in " << dim << "D with degree " << fe_degree << "..."
          << std::endl;
  // 1. Generate triangulation
  Triangulation<dim> triangulation(
    Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(); // Refine once

  // 2. Initialize FE, DoFHandler, Mapping

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  const unsigned int n_lanes = VectorizedArray<double>::size();

  MappingQ1<dim> mapping; // Using Q1 mapping as in step-94

  // 3. Initialize constraints (needed for MatrixFree)
  AffineConstraints<double> constraints;
  IndexSet                  locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
  constraints.reinit(locally_relevant_dofs);
  // No boundary conditions or hanging nodes needed for this simple test yet
  constraints.close();

  // 4. Initialize MatrixFree for level 1
  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mg_level          = 1; // Target level 1
  additional_data.store_ghost_cells = true;
  // Add mapping flags needed for FEEvaluation/FEPatchEvaluation
  additional_data.mapping_update_flags =
    (update_gradients | update_JxW_values | update_quadrature_points);

  std::shared_ptr<MatrixFree<dim, double>> mf_level_storage =
    std::make_shared<MatrixFree<dim, double>>();
  mf_level_storage->reinit(mapping,
                           dof_handler,
                           constraints,
                           QGauss<1>(fe_degree + 1),
                           additional_data);

  // 5. Build PatchStorage on level 1
  PatchStorage<MatrixFree<dim, double>> patch_storage(mf_level_storage);
  patch_storage.initialize();

  // 6. Initialize vectors
  LinearAlgebra::distributed::Vector<double> src, dst_mf, dst_patch;
  mf_level_storage->initialize_dof_vector(src);
  mf_level_storage->initialize_dof_vector(dst_mf);
  mf_level_storage->initialize_dof_vector(dst_patch);

  // Fill src with random values
  for (unsigned int i = 0; i < src.locally_owned_size(); ++i)
    src.local_element(i) = i;
  src.update_ghost_values();

  // 7. Apply Laplace via standard MatrixFree

  using FEEval = FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double>;
  FEEval phi(*mf_level_storage);
  for (unsigned int cell = 0; cell < mf_level_storage->n_cell_batches(); ++cell)
    {
      phi.reinit(cell);
      phi.read_dof_values(src);
      phi.evaluate(EvaluationFlags::gradients);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_gradient(phi.get_gradient(q), q);
      phi.integrate(EvaluationFlags::gradients);
      phi.distribute_local_to_global(dst_mf);
    }


  dst_mf.compress(VectorOperation::add);


  // 8. Apply Laplace via FEPatchEvaluation

  using PatchEval =
    FEPatchEvaluation<FEEval, PatchDistributors::Dynamic<dim, fe_degree>>;

  // Test constructor and  move constructor: consistency of cell_dofs_view_raw
  // can fail. If so, it will be detected in reinit(). PatchEval
  // patch_eval(patch_storage, FEEval(*mf_level_storage));
  PatchEval patch_eval(
    [&]() { return PatchEval(patch_storage, FEEval(*mf_level_storage)); }());

  // Checks: there should only one patch

  Assert(patch_storage.n_patches() > 0, ExcInternalError());
  deallog << "  Level 1 PatchStorage initialized with "
          << patch_storage.n_patches() << " patches." << std::endl;

  patch_eval.reinit(0);
  patch_eval.read_dof_values(src);

  for (auto &phi : patch_eval.fe_evaluations)
    {
      phi.evaluate(EvaluationFlags::gradients);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_gradient(phi.get_gradient(q), q);
      phi.integrate(EvaluationFlags::gradients);
    }

  // Allocate patch vectors
  AlignedVector<double> patch_result_local(patch_eval.n_patch_dofs());
  AlignedVector<double> patch_result_global(
    patch_eval.n_patch_dofs()); // Second vector for later use if needed

  // Gather results from FEEvaluation objects into the first patch vector
  patch_eval.gather_local_to_patch(ArrayView<double>(patch_result_local),
                                   true); // Accumulate results

  // Distribute results from cell evaluations to the global vector dst_patch

  patch_eval.read_dof_values(dst_mf);
  patch_eval.gather_local_to_patch(ArrayView<double>(patch_result_global),
                                   false);

  // 9. Compare patch_result_local (accumulated) and patch_result_global
  // (non-accumulated)
  AssertDimension(patch_result_local.size(), patch_result_global.size());

  for (unsigned int i = 0; i < patch_result_local.size(); ++i)
    {
      const double diff = patch_result_local[i] - patch_result_global[i];

      if (diff * diff > 1e-10)
        deallog << "Difference " << diff << " :" << patch_result_local[i]
                << "  " << patch_result_global[i] << " on DoF " << i
                << std::endl;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  initlog();

  test<2, 1>();
  test<2, 2>();
  test<2, 3>();

  test<2, 4>();
  test<2, 5>();
  test<2, 6>();
  test<2, 7>();


  test<3, 1>();
  test<3, 2>();
  test<3, 3>();

  test<3, 4>();
  test<3, 5>();
  test<3, 6>();
  test<3, 7>();
  deallog << "Tests finished." << std::endl;

  return 0; // Indicate success
}
