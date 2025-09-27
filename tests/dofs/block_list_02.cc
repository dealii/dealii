// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include "block_list.h"

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
    deallog.push("ttttt");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, true, true, true, true, true);
    bl.compress();
    print_patches(bl);
    deallog.push("inv vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("ttfft");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, true, true, false, false, true);
    bl.compress();
    print_patches(bl);
    deallog.push("inv vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("tffft");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, true, false, false, false, true);
    bl.compress();
    print_patches(bl);
    deallog.push("inv vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("ftttt");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, false, true, true, true, true);
    bl.compress();
    print_patches(bl);
    deallog.push("inv vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("fffft");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, false, false, false, false, true);
    bl.compress();
    print_patches(bl);
    deallog.push("inv vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("ttttf");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, true, true, true, true, false);
    deallog.push("vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("ttfff");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, true, true, false, false, false);
    deallog.push("vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("tffff");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, true, false, false, false, false);
    deallog.push("vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("ftttf");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, false, true, true, true, false);
    deallog.push("vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("fffff");
    SparsityPattern           bl;
    std::vector<unsigned int> vm;
    vm = DoFTools::make_vertex_patches(
      bl, dof, level, false, false, false, false, false);
    deallog.push("vertex mapping");
    for (unsigned int i = 0; i < vm.size(); ++i)
      deallog << ' ' << vm[i];
    deallog << std::endl;
    deallog.pop();
    deallog.pop();
    deallog << std::endl;
  }
}


int
main()
{
  initlog();
  deallog.push("2D");
  test_global_refinement<Triangulation<2>>(&test_block_list<2>);
  deallog.pop();
  deallog.push("3D");
  test_global_refinement<Triangulation<3>>(&test_block_list<3>);
  deallog.pop();
}
