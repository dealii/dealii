// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// this is a template for matrix-vector products with the Helmholtz equation
// (zero and first derivatives) on different kinds of meshes (Cartesian,
// general, with and without hanging nodes). It also tests the multithreading

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


// forward declare this function. will be implemented in .cc files
template <int dim, int fe_degree, typename Number>
void
test();

template <int dim>
double
ConstantCoefficient(const Point<dim> & /*p*/)
{
  return 10.;
}


template <int dim>
double
VaryingCoefficient(const Point<dim> &p)
{
  return 10. / (0.05 + 2. * p.square());
}


template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType,
          int n_q_points_1d>
void
do_test(const DoFHandler<dim>           &dof,
        const AffineConstraints<Number> &constraints,
        const unsigned int               n_locally_owned_cells,
        const bool                       constant_coefficient = true,
        const bool                       coloring             = false)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MappingQ<dim>                                              mapping(fe_degree);
  Portable::MatrixFree<dim, Number>                          mf_data;
  typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients |
                                         update_JxW_values |
                                         update_quadrature_points;
  additional_data.use_coloring = coloring;
  const QGauss<1> quad(n_q_points_1d);
  mf_data.reinit(mapping, dof, constraints, quad, additional_data);

  const unsigned int n_dofs = dof.n_dofs();
  MatrixFreeTest<dim, fe_degree, Number, VectorType, n_q_points_1d> mf(
    mf_data,
    n_locally_owned_cells * std::pow(n_q_points_1d, dim),
    constant_coefficient);
  Vector<Number>                         in_host(n_dofs), out_host(n_dofs);
  LinearAlgebra::ReadWriteVector<Number> in(n_dofs), out(n_dofs);
  VectorType                             in_device(n_dofs);
  VectorType                             out_device(n_dofs);

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
    QGauss<dim> quadrature_formula(n_q_points_1d);

    FEValues<dim> fe_values(mapping,
                            dof.get_fe(),
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                   endc = dof.end();
    for (; cell != endc; ++cell)
      {
        cell_matrix = 0;
        fe_values.reinit(cell);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const auto coef =
              constant_coefficient ?
                ConstantCoefficient<dim>(fe_values.quadrature_point(q_point)) :
                VaryingCoefficient<dim>(fe_values.quadrature_point(q_point));
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

  Number out_norm = 0.;
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
    deallog.push("double");
    deallog.push("2d");
    test<2, 1, double>();
    test<2, 2, double>();
    test<2, 3, double>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1, double>();
    test<3, 2, double>();
    test<3, 3, double>();
    deallog.pop();
    deallog.pop();
    deallog.push("float");
    deallog.push("2d");
    test<2, 1, double>();
    test<2, 2, double>();
    test<2, 3, double>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1, double>();
    test<3, 2, double>();
    test<3, 3, double>();
    deallog.pop();
    deallog.pop();
  }

  Kokkos::finalize();

  return 0;
}
