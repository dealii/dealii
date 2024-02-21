// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

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

template <int dim>
void
test_block_list(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs();

  const unsigned int level = tr.n_levels() - 1;

  {
    deallog.push("(tt)ffff");
    SparsityPattern           bl;
    std::vector<bool>         temp_vector(2, true);
    BlockMask                 exclude_boundary_dofs(temp_vector);
    std::vector<unsigned int> vm;
    std::cout << exclude_boundary_dofs.size() << std::endl;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, exclude_boundary_dofs, false, false, false, false);
    bl.compress();
    print_patches(bl);
    deallog.push("vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("(tf)ffff");
    SparsityPattern   bl;
    std::vector<bool> temp_vector(2, true);
    temp_vector[1] = false;
    BlockMask                 exclude_boundary_dofs(temp_vector);
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, exclude_boundary_dofs, false, false, false, false);
    bl.compress();
    print_patches(bl);
    deallog.push("vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("(ft)ffff");
    SparsityPattern   bl;
    std::vector<bool> temp_vector(2, false);
    temp_vector[1] = true;
    BlockMask                 exclude_boundary_dofs(temp_vector);
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, exclude_boundary_dofs, false, false, false, false);
    bl.compress();
    print_patches(bl);
    deallog.push("vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("(ff)ffff");
    SparsityPattern           bl;
    std::vector<bool>         temp_vector(2, false);
    BlockMask                 exclude_boundary_dofs(temp_vector);
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, exclude_boundary_dofs, false, false, false, false);
    bl.compress();
    print_patches(bl);
    deallog.push("vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
}

template <int dim>
void
test_global_refinement(void (*test_block_list)(const Triangulation<dim> &tr,
                                               const FiniteElement<dim> &fe))
{
  Triangulation<dim> trc(
    Triangulation<dim>::limit_level_difference_at_vertices);
  Triangulation<dim> trl(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(trc);
  trc.refine_global(2);
  GridGenerator::hyper_L(trl);
  trl.refine_global(1);

  FESystem<dim, dim> fe1(FESystem<dim, dim>(FE_Q<dim, dim>(2), dim),
                         1,
                         FE_Q<dim, dim>(1),
                         1);

  deallog.push("Square");
  test_block_list(trc, fe1);
  deallog.pop();
  deallog.push("L");
  test_block_list(trl, fe1);
  deallog.pop();
}



int
main()
{
  initlog();
  deallog.push("2D");
  test_global_refinement<2>(&test_block_list<2>);
  deallog.pop();
  deallog.push("3D");
  test_global_refinement<3>(&test_block_list<3>);
  deallog.pop();
}
