// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_transfer.h>

#include <string>
#include <vector>

#include "../tests.h"


#define TEST(dim, l, el, deg)                                        \
  {                                                                  \
    el<dim> fe(deg);                                                 \
    deallog << #el << '<' << dim << ">(" << deg << ')' << std::endl; \
    print_matrix(tr##dim, l, fe, #el);                               \
  }

template <int dim>
inline void
print_matrix(Triangulation<dim>       &tr,
             unsigned int              level,
             const FiniteElement<dim> &finel,
             const char * /*name*/)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(finel);
  dof.distribute_mg_dofs();

  MGTransferPrebuilt<Vector<double>> transfer;
  transfer.build(dof);

  unsigned int   n_coarse = dof.n_dofs(level - 1);
  unsigned int   n_fine   = dof.n_dofs(level);
  Vector<double> in(n_fine);
  Vector<double> out(n_coarse);

  for (unsigned int i = 0; i < n_fine; ++i)
    {
      in    = 0.;
      out   = 0.;
      in(i) = 1.;
      transfer.restrict_and_add(level, out, in);
      for (unsigned int k = 0; k < out.size(); ++k)
        deallog << '\t' << out(k);
      deallog << std::endl;
    }
  deallog << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(7);

  Triangulation<2> tr2(Triangulation<2>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(tr2, -1., 1.);
  tr2.refine_global(2);

  Triangulation<3> tr3(Triangulation<3>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(tr3, -1., 1.);
  tr3.refine_global(3);

  TEST(2, 1, FE_Q, 1);
  TEST(2, 1, FE_Q, 2);
  TEST(2, 1, FE_Q, 3);
  //  TEST(2, 1, FE_Q, 4);

  TEST(2, 1, FE_DGQ, 0);
  TEST(2, 1, FE_DGQ, 1);
  TEST(2, 1, FE_DGQ, 2);
  TEST(2, 1, FE_DGQ, 3);
  TEST(2, 1, FE_DGQ, 4);

  TEST(2, 1, FE_DGP, 0);
  TEST(2, 1, FE_DGP, 1);
  TEST(2, 1, FE_DGP, 2);
  TEST(2, 1, FE_DGP, 3);
  TEST(2, 1, FE_DGP, 4);
  TEST(2, 1, FE_DGP, 5);
  TEST(2, 1, FE_DGP, 6);

  TEST(3, 1, FE_DGP, 0);
  TEST(3, 1, FE_DGP, 1);
  TEST(3, 1, FE_DGP, 2);
  TEST(3, 1, FE_DGP, 3);
  TEST(3, 1, FE_DGP, 4);

  return 0;
}
