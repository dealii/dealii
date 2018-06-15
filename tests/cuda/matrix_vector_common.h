// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// this is a template for matrix-vector products with the Helmholtz equation
// (zero and first derivatives) on different kinds of meshes (Cartesian,
// general, with and without hanging nodes). It also tests the multithreading

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"
#include "matrix_vector_mf.h"


// forward declare this function. will be implemented in .cc files
template <int dim, int fe_degree>
void
test();



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
void
do_test(const DoFHandler<dim> &          dof,
        const AffineConstraints<double> &constraints)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MappingQGeneric<dim>                  mapping(fe_degree);
  CUDAWrappers::MatrixFree<dim, Number> mf_data;
  typename CUDAWrappers::MatrixFree<dim, Number>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients |
                                         update_JxW_values |
                                         update_quadrature_points;
  const QGauss<1> quad(n_q_points_1d);
  mf_data.reinit(mapping, dof, constraints, quad, additional_data);

  const unsigned int                                    n_dofs = dof.n_dofs();
  MatrixFreeTest<dim, fe_degree, Number, n_q_points_1d> mf(mf_data);
  Vector<Number>                              in_host(n_dofs), out_host(n_dofs);
  LinearAlgebra::ReadWriteVector<Number>      in(n_dofs), out(n_dofs);
  LinearAlgebra::CUDAWrappers::Vector<Number> in_device(n_dofs);
  LinearAlgebra::CUDAWrappers::Vector<Number> out_device(n_dofs);

  for (unsigned int i = 0; i < n_dofs; ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = Testing::rand() / (double)RAND_MAX;
      in(i)              = entry;
    }

  in_device.import(in, VectorOperation::insert);
  mf.vmult(out_device, in_device);
  cudaDeviceSynchronize();
  out.import(out_device, VectorOperation::insert);


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
                              update_JxW_values);

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
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) += ((fe_values.shape_grad(i, q_point) *
                                         fe_values.shape_grad(j, q_point) +
                                       10. * fe_values.shape_value(i, q_point) *
                                         fe_values.shape_value(j, q_point)) *
                                      fe_values.JxW(q_point));
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

  Number out_dist_cpu_norm = 0.;
  Number out_norm          = 0.;
  for (unsigned i = 0; i < n_dofs; ++i)
    {
      out_norm += std::pow(out[i] - out_host[i], 2);
      out_dist_cpu_norm += std::pow(out_host[i], 2);
    }
  const double diff_norm = out_norm / out_dist_cpu_norm;
  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}



int
main()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << std::setprecision(3);

  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    test<2, 3>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    test<3, 3>();
    deallog.pop();
  }

  return 0;
}
