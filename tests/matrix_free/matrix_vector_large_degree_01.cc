// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that FEEvaluation can correctly handle large polynomial degrees of 20
// in 3D. The old implementation of FEEvaluation that allocated all data
// needed by the kernel on the stack ran into stack overflows in that case
// since degree 20 needs around 3.5 MB of data. Since we cannot compare with
// matrix-based results in this case, this test simply checks that the
// matrix-vector product runs without error

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

  Vector<double> in(dof.n_dofs()), out(dof.n_dofs());

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      in(i) = 1.;
    }

  {
    MatrixFree<dim, double> mf_data;
    {
      const QGauss<1>                                  quad(fe_degree + 1);
      typename MatrixFree<dim, double>::AdditionalData data;
      data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::none;
      mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
    }

    MatrixFreeTest<dim, fe_degree, double, Vector<double>, fe_degree + 1> mf(
      mf_data);
    mf.vmult(out, in);
    deallog << "Result norm: " << out.l2_norm() << std::endl;
  }
  {
    MatrixFree<dim, double> mf_data;
    {
      const QGauss<1>                                  quad(fe_degree + 2);
      typename MatrixFree<dim, double>::AdditionalData data;
      data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::none;
      mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
    }

    MatrixFreeTest<dim, fe_degree, double, Vector<double>, fe_degree + 2> mf(
      mf_data);
    mf.vmult(out, in);
    deallog << "Result norm: " << out.l2_norm() << std::endl;
  }
}



int
main()
{
  initlog();

  {
    deallog.push("2d");
    test<2, 5>();
    test<2, 15>();
    test<2, 35>();
    deallog.pop();
    deallog.push("3d");
    test<3, 10>();
    test<3, 20>();
    deallog.pop();
  }
}
