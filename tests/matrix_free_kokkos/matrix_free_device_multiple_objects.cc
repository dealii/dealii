// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Reproduce a bug where only one matrix-free object is valid

#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

#include "Kokkos_Core.hpp"
#include "matrix_vector_device_mf.h"

template <int fe_degree, int n_q_points_1d>
void
do_test(const DoFHandler<2> &dof,
        MatrixFreeTest<
          2,
          fe_degree,
          double,
          LinearAlgebra::distributed::Vector<double, MemorySpace::Default>,
          n_q_points_1d>                &mf,
        unsigned int                     n_dofs,
        MappingQ<2>                     &mapping,
        const AffineConstraints<double> &constraints)
{
  Vector<double>                         in_host(n_dofs), out_host(n_dofs);
  LinearAlgebra::ReadWriteVector<double> in(n_dofs), out(n_dofs);
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> in_device(
    n_dofs);
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> out_device(
    n_dofs);

  for (unsigned int i = 0; i < n_dofs; ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = Testing::rand() / (double)RAND_MAX;
      in(i)              = entry;
    }

  in_device.import_elements(in, VectorOperation::insert);
  mf.vmult(out_device, in_device);
  Kokkos::fence();
  out.import_elements(out_device, VectorOperation::insert);

  // assemble sparse matrix with (\nabla v, \nabla u) + (v, 10 * u)
  SparsityPattern sparsity;
  {
    DynamicSparsityPattern csp(n_dofs, n_dofs);
    DoFTools::make_sparsity_pattern(dof, csp, constraints, true);
    sparsity.copy_from(csp);
  }
  SparseMatrix<double> sparse_matrix(sparsity);
  {
    QGauss<2> quadrature_formula(n_q_points_1d);

    FEValues<2> fe_values(mapping,
                          dof.get_fe(),
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<2>::active_cell_iterator cell = dof.begin_active(),
                                                 endc = dof.end();
    for (; cell != endc; ++cell)
      {
        cell_matrix = 0;
        fe_values.reinit(cell);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const auto coef = 10.;
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    ((fe_values.shape_grad(i, q_point) *
                        fe_values.shape_grad(j, q_point) +
                      coef * fe_values.shape_value(i, q_point) *
                        fe_values.shape_value(j, q_point)) *
                     fe_values.JxW(q_point));
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               local_dof_indices,
                                               sparse_matrix);
      }
  }
  for (unsigned i = 0; i < n_dofs; ++i)
    in_host[i] = in[i];
  sparse_matrix.vmult(out_host, in_host);

  double out_norm = 0.;
  for (unsigned i = 0; i < n_dofs; ++i)
    out_norm += std::pow(out[i] - out_host[i], 2);
  const double diff_norm = std::sqrt(out_norm) / out_host.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}


int
main()
{
  initlog();
  deallog.depth_console(0);

  deallog << std::setprecision(3);

  Kokkos::initialize();

  {
    Triangulation<2> tria;
    GridGenerator::hyper_cube(tria);
    tria.refine_global(5 - 2);
    AffineConstraints<double> constraints;
    constraints.close();
    bool constant_coefficient = true;

    // Create the first MatrixFree object
    constexpr unsigned int fe_degree_1     = 1;
    constexpr unsigned int n_q_points_1d_1 = fe_degree_1 + 1;
    FE_Q<2>                fe_1(fe_degree_1);
    DoFHandler<2>          dof_1(tria);
    dof_1.distribute_dofs(fe_1);
    MappingQ<2>                                     mapping_1(fe_degree_1);
    Portable::MatrixFree<2, double>                 mf_data_1;
    Portable::MatrixFree<2, double>::AdditionalData additional_data_1;
    additional_data_1.mapping_update_flags = update_values | update_gradients |
                                             update_JxW_values |
                                             update_quadrature_points;
    const QGauss<1> quad_1(n_q_points_1d_1);
    mf_data_1.reinit(mapping_1, dof_1, constraints, quad_1, additional_data_1);
    const unsigned int n_dofs_1 = dof_1.n_dofs();
    MatrixFreeTest<
      2,
      fe_degree_1,
      double,
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default>,
      n_q_points_1d_1>
      mf_1(mf_data_1,
           n_dofs_1 * std::pow(n_q_points_1d_1, 2),
           constant_coefficient);

    // Create the second MatrixFree object
    constexpr unsigned int fe_degree_2     = 2;
    constexpr unsigned int n_q_points_1d_2 = fe_degree_2 + 1;
    FE_Q<2>                fe_2(fe_degree_2);
    DoFHandler<2>          dof_2(tria);
    dof_2.distribute_dofs(fe_2);
    MappingQ<2>                                     mapping_2(fe_degree_2);
    Portable::MatrixFree<2, double>                 mf_data_2;
    Portable::MatrixFree<2, double>::AdditionalData additional_data_2;
    additional_data_2.mapping_update_flags = update_values | update_gradients |
                                             update_JxW_values |
                                             update_quadrature_points;
    const QGauss<1> quad_2(n_q_points_1d_2);
    mf_data_2.reinit(mapping_2, dof_2, constraints, quad_2, additional_data_2);
    const unsigned int n_dofs_2 = dof_2.n_dofs();
    MatrixFreeTest<
      2,
      fe_degree_2,
      double,
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default>,
      n_q_points_1d_2>
      mf_2(mf_data_2,
           n_dofs_2 * std::pow(n_q_points_1d_2, 2),
           constant_coefficient);

    // Perform MV with the first object
    do_test<fe_degree_1, n_q_points_1d_1>(
      dof_1, mf_1, n_dofs_1, mapping_1, constraints);

    // Perform MV with the second object
    do_test<fe_degree_2, n_q_points_1d_2>(
      dof_2, mf_2, n_dofs_2, mapping_2, constraints);
  }
  Kokkos::finalize();

  return 0;
}
