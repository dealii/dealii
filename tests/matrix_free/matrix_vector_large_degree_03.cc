// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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



// test that FEEvaluation can correctly handle large polynomial degrees of 20
// with 100, 200, 500 and 1000 quadrature points per direction in 2D. Since we
// have a constant-coefficient Laplacian, we can verify the implementation by
// comparing to the result with 21 quadrature points.

#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "matrix_vector_mf.h"


template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  Vector<double> in(dof.n_dofs()), out(dof.n_dofs()), ref(dof.n_dofs());

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      in(i) = static_cast<double>(i % 13) / 11.;
    }

  {
    MatrixFree<dim, double> mf_data;
    {
      const QGauss<1>                                  quad(fe_degree + 1);
      typename MatrixFree<dim, double>::AdditionalData data;
      data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::none;
      mf_data.reinit(dof, constraints, quad, data);
    }

    MatrixFreeTest<dim, fe_degree, double, Vector<double>, fe_degree + 1> mf(
      mf_data);
    mf.vmult(ref, in);
  }

  typename MatrixFree<dim, double>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
  MatrixFree<dim, double> mf_data;

  {
    mf_data.reinit(dof, constraints, QGauss<1>(fe_degree + 2), data);
    MatrixFreeTest<dim, fe_degree, double, Vector<double>, fe_degree + 2> mf(
      mf_data);
    mf.vmult(out, in);
    out -= ref;
    deallog << "Error with " << fe_degree + 2 << "^" << dim
            << " quadrature points: " << out.l2_norm() << std::endl;
  }

  // unfortunately we cannot use for loops due to the template, so duplicate
  // some code here
  {
    mf_data.reinit(dof, constraints, QGauss<1>(100), data);
    MatrixFreeTest<dim, fe_degree, double, Vector<double>, 100> mf(mf_data);
    mf.vmult(out, in);
    out -= ref;
    deallog << "Error with " << 100 << "^" << dim
            << " quadrature points: " << out.l2_norm() << std::endl;
  }
  {
    mf_data.reinit(dof, constraints, QGauss<1>(200), data);
    MatrixFreeTest<dim, fe_degree, double, Vector<double>, 200> mf(mf_data);
    mf.vmult(out, in);
    out -= ref;
    deallog << "Error with " << 200 << "^" << dim
            << " quadrature points: " << out.l2_norm() << std::endl;
  }
  {
    mf_data.reinit(dof, constraints, QGauss<1>(500), data);
    MatrixFreeTest<dim, fe_degree, double, Vector<double>, 500> mf(mf_data);
    mf.vmult(out, in);
    out -= ref;
    deallog << "Error with " << 500 << "^" << dim
            << " quadrature points: " << out.l2_norm() << std::endl;
  }
  {
    mf_data.reinit(dof, constraints, QGauss<1>(1000), data);
    MatrixFreeTest<dim, fe_degree, double, Vector<double>, 1000> mf(mf_data);
    mf.vmult(out, in);
    out -= ref;
    deallog << "Error with " << 1000 << "^" << dim
            << " quadrature points: " << out.l2_norm() << std::endl;
  }
}



int
main()
{
  initlog();

  test<2, 20>();
}
