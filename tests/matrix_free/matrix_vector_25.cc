// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// similar test as matrix_vector_01, but using ArrayView instead of
// dealii::Vector

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
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

#include "matrix_vector_mf.h"


template <int dim, int fe_degree>
void
test();



template <int dim, typename number>
void
do_test(const DoFHandler<dim>           &dof,
        const AffineConstraints<double> &constraints)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MappingQ<dim>           mapping(dof.get_fe().degree);
  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1> quad(dof.get_fe().degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::partition_partition;

    mf_data.reinit(mapping, dof, constraints, quad, data);
  }

  Vector<number> in(dof.n_dofs()), out(dof.n_dofs());
  Vector<number> in_dist(dof.n_dofs());
  Vector<number> out_dist(in_dist);

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = random_value<double>();
      in(i)              = entry;
      in_dist(i)         = entry;
    }

  MatrixFreeTest<dim, -1, number, Vector<number>, 0> mf_reference(mf_data);
  mf_reference.vmult(out, in);

  MatrixFreeTest<dim, -1, number, dealii::ArrayView<number>, 0> mf(mf_data);
  ArrayView<number> out_view = make_array_view(out_dist);
  mf.vmult(out_view, make_array_view(in_dist));


  {
    out_dist -= out;
    const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
    deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
  }

  mf.vmult(out_view, make_array_view(in_dist));

  {
    out_dist -= out;
    const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
    deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
  }
}

template <int dim>
void
test(const unsigned int fe_degree)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5 - dim);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  do_test<dim, double>(dof, constraints);
}



int
main()
{
  initlog();

  {
    deallog.push("2d");
    test<2>(1);
    test<2>(2);
    deallog.pop();
    deallog.push("3d");
    test<3>(2);
    deallog.pop();
  }
}
