// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <algorithm>

#include "../tests.h"



void
print_patches(const SparsityPattern &bl)
{
  for (unsigned int i = 0; i < bl.n_rows(); ++i)
    {
      deallog << "Block " << std::setw(3) << i;
      std::vector<unsigned int> entries;
      for (SparsityPattern::iterator b = bl.begin(i); b != bl.end(i); ++b)
        entries.push_back(b->column());

      std::sort(entries.begin(), entries.end());

      for (unsigned int i = 0; i < entries.size(); ++i)
        deallog << ' ' << std::setw(4) << entries[i];
      deallog << std::endl;
    }
}


template <class TR>
void
test_global_refinement(
  void (*test_block_list)(const TR &tr, const FiniteElement<TR::dimension> &fe))
{
  const unsigned int dim = TR::dimension;
  TR trc(Triangulation<dim>::limit_level_difference_at_vertices);
  TR trl(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(trc);
  trc.refine_global(2);
  GridGenerator::hyper_L(trl);
  trl.refine_global(1);

  FE_DGQ<dim> fe1(0);
  FE_DGQ<dim> fe2(1);
  //  FE_Q<dim> fe3(1);
  FE_RaviartThomas<dim> fe4(0);
  FE_RaviartThomas<dim> fe5(1);
  FE_Nedelec<dim>       fe6(0);
  FE_Nedelec<dim>       fe7(1);

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

template <int dim>
void
test_global_refinement_parallel(
  void (*test_block_list)(const parallel::distributed::Triangulation<dim> &tr,
                          const FiniteElement<dim>                        &fe))
{
  parallel::distributed::Triangulation<dim> trl(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_L(trl);
  trl.refine_global(2);

  FE_DGQ<dim>           fe1(0);
  FE_DGQ<dim>           fe2(1);
  FE_Q<dim>             fe3(1);
  FE_RaviartThomas<dim> fe4(0);
  FE_Nedelec<dim>       fe6(0);

  deallog.push("L");
  test_block_list(trl, fe1);
  test_block_list(trl, fe2);
  test_block_list(trl, fe3);
  test_block_list(trl, fe4);
  test_block_list(trl, fe6);
  deallog.pop();
}
