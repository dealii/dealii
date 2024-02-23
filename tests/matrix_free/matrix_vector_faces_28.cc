// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that MatrixFree::update_mapping() works as expected. The test is
// otherwise similar to matrix_vector_faces_03.

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>

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
  // if (fe_degree == 1)
  //  tria.refine_global(1);
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

  deallog << "Degree of element: " << fe_degree << std::endl;

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  FESystem<dim>   fe_grid(FE_Q<dim>(dof.get_fe().degree), dim);
  DoFHandler<dim> dof_grid(tria);
  dof_grid.distribute_dofs(fe_grid);
  Vector<double> euler;
  euler.reinit(dof_grid.n_dofs());
  const ComponentMask mask(dim, true);
  VectorTools::get_position_vector(dof_grid, euler, mask);

  MappingFEField<dim> mapping(dof_grid, euler, mask);

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

  auto create_reference_result = [&]() {
    // assemble sparse matrix with MeshWorker
    SparsityPattern      sparsity;
    SparseMatrix<double> matrix;
    {
      DynamicSparsityPattern d_sparsity(dof.n_dofs());
      DoFTools::make_flux_sparsity_pattern(dof, d_sparsity);
      sparsity.copy_from(d_sparsity);
    }
    const unsigned int n_q_points_1d = dof.get_fe().degree + 1;
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
  };

  MatrixFree<dim>                          mf_data;
  const QGauss<1>                          quad(dof.get_fe().degree + 1);
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  data.tasks_block_size      = 3;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  mf_data.reinit(mapping, dof, constraints, quad, data);

  MatrixFreeTest<dim, fe_degree, fe_degree + 1, double, Vector<double>, 1> mf(
    mf_data);
  mf.vmult(out_dist, in);
  create_reference_result();

  out_dist -= out;
  double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference:           " << diff_norm << std::endl;

  // deform the geometry slightly
  euler(0) += 0.001;
  euler(euler.size() - 1) += 0.001;

  mf.vmult(out_dist, in);
  create_reference_result();

  out_dist -= out;
  diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference no update: " << diff_norm << std::endl;

  mf_data.update_mapping(mapping);
  mf.vmult(out_dist, in);
  out_dist -= out;
  diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference w update:  " << diff_norm << std::endl;
}
