// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_nedelec.h>

#include <algorithm>

using namespace dealii;


void
print_patches (const SparsityPattern &bl)
{
  for (unsigned int i=0; i<bl.n_rows(); ++i)
    {
      deallog << "Block " << std::setw(3) << i;
      std::vector<unsigned int> entries;
      for (SparsityPattern::row_iterator b = bl.row_begin(i); b != bl.row_end(i); ++b)
        entries.push_back(*b);

      std::sort(entries.begin(), entries.end());

      for (unsigned int i=0; i<entries.size(); ++i)
        deallog << ' ' << std::setw(4) << entries[i];
      deallog << std::endl;
    }
}


template <int dim>
void test_global_refinement(
  void (*test_block_list)(const Triangulation<dim> &tr, const FiniteElement<dim> &fe))
{
  Triangulation<dim> trc, trl;
  GridGenerator::hyper_cube(trc);
  trc.refine_global(2);
  GridGenerator::hyper_L(trl);
  trl.refine_global(1);

  FE_DGQ<dim> fe1(0);
  FE_DGQ<dim> fe2(1);
//  FE_Q<dim> fe3(1);
  FE_RaviartThomas<dim> fe4(0);
  FE_RaviartThomas<dim> fe5(1);
  FE_Nedelec<dim> fe6(0);
  FE_Nedelec<dim> fe7(1);

  deallog.push("Square");
  test_block_list(trc, fe1);
  test_block_list(trc, fe2);
//  test_block_list(trc, fe3);
  test_block_list(trc, fe4);
  test_block_list(trc, fe5);
  test_block_list(trc, fe6);
  test_block_list(trc, fe7);
  deallog.pop();
  deallog.push("L");
  test_block_list(trl, fe1);
  test_block_list(trl, fe2);
//  test_block_list(trl, fe3);
  test_block_list(trl, fe4);
  test_block_list(trl, fe5);
  test_block_list(trl, fe6);
  test_block_list(trl, fe7);
  deallog.pop();
}
