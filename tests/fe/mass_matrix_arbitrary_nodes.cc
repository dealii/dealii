// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// verify that we get a diagonal matrix when using FE_DGQArbitraryNodes with
// the same quadrature formula for integration as for construction of the FE.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/numerics/matrix_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim>
void
check()
{
  deallog << dim << 'D' << std::endl;
  for (unsigned int q_points = 1; q_points < 7 - dim; ++q_points)
    {
      deallog << "q_points=" << q_points << std::endl;

      Triangulation<dim> tr;
      if (dim == 1)
        GridGenerator::hyper_cube(tr);
      else
        GridGenerator::hyper_ball(tr);
      tr.reset_manifold(0);
      tr.refine_global(3 - dim);

      FE_DGQArbitraryNodes<dim> fe((QGauss<1>(q_points)));
      DoFHandler<dim>           dh(tr);
      dh.distribute_dofs(fe);

      SparsityPattern sp(dh.n_dofs(), dh.n_dofs(), fe.dofs_per_cell);
      DoFTools::make_sparsity_pattern(dh, sp);
      sp.compress();

      SparseMatrix<double> mass_matrix(sp);
      MatrixTools::create_mass_matrix(dh, QGauss<dim>(q_points), mass_matrix);

      // verify that the matrix is diagonal
      for (unsigned int i = 0; i < dh.n_dofs(); ++i)
        {
          deallog << mass_matrix.el(i, i) << std::endl;
          for (unsigned int j = 0; j < dh.n_dofs(); ++j)
            if (i != j)
              Assert(std::fabs(mass_matrix.el(i, j) / mass_matrix.el(i, i)) <
                       1e-14,
                     ExcInternalError());
        }
    }
}



int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
}
