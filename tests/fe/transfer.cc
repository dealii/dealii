//----------------------------  show_transfer.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  show_transfer.cc  ---------------------------
//
// Print multigrid transfer matrices between one and four cells.
//
//----------------------------  show_transfer.cc  ---------------------------

#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_values.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_dof_handler.h>

#include <vector>
#include <fstream>
#include <string>

char fname[50];

#define TEST(dim, l, el, deg) { el<dim> fe(deg); \
  of << # el << '<' << dim << ">(" << deg << ')' << std::endl; \
  print_matrix(of, tr ## dim, l, fe, #el); }

template<int dim>
inline void
print_matrix(std::ostream& of,
	     Triangulation<dim>& tr,
	     unsigned int level,
	     const FiniteElement<dim>& finel,
	     const char* name)
{
  MGDoFHandler<dim> dof(tr);
  dof.distribute_dofs(finel);

  MGTransferPrebuilt<double> transfer;
  transfer.build_matrices(dof);

  unsigned int n_coarse = dof.n_dofs(level-1);
  unsigned int n_fine = dof.n_dofs(level);
  Vector<double> in(n_fine);
  Vector<double> out(n_coarse);

  for (unsigned int i=0;i<n_fine;++i)
    {
      in = 0.;
      out = 0.;
      in(i) = 1.;
      transfer.restrict_and_add(level,out,in);
      out.print(of, 3, false);
    }
  of << std::endl;
}


int
main()
{
  Triangulation<2> tr2;

  GridGenerator::hyper_cube(tr2, -1., 1.);
  tr2.refine_global(2);

  Triangulation<3> tr3;

  GridGenerator::hyper_cube(tr3, -1., 1.);
  tr3.refine_global(3);

  std::ofstream of("transfer.output");
  
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

