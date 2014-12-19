// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <vector>
#include <fstream>
#include <string>


#define TEST(dim, l, el, deg) { el<dim> fe(deg); \
    deallog << # el << '<' << dim << ">(" << deg << ')' << std::endl; \
    print_matrix(tr ## dim, l, fe, #el); }

template<int dim>
inline void
print_matrix(Triangulation<dim> &tr,
             unsigned int level,
             const FiniteElement<dim> &finel,
             const char * /*name*/)
{
  MGDoFHandler<dim> dof(tr);
  dof.distribute_dofs(finel);

  MGTransferPrebuilt<Vector<double> > transfer;
  transfer.build_matrices(dof);

  unsigned int n_coarse = dof.n_dofs(level-1);
  unsigned int n_fine = dof.n_dofs(level);
  Vector<double> in(n_fine);
  Vector<double> out(n_coarse);

  for (unsigned int i=0; i<n_fine; ++i)
    {
      in = 0.;
      out = 0.;
      in(i) = 1.;
      transfer.restrict_and_add(level,out,in);
      for (unsigned int k=0; k<out.size(); ++k)
        deallog << '\t' << out(k);
      deallog << std::endl;
    }
  deallog << std::endl;
}


int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision(7);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> tr2;

  GridGenerator::hyper_cube(tr2, -1., 1.);
  tr2.refine_global(2);

  Triangulation<3> tr3;

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
