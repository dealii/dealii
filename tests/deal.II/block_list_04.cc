// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
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


#include "block_list.h"

template <int dim>
void
test_block_list(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  DoFHandler<dim> dof;
  dof.initialize(tr, fe);
  dof.distribute_mg_dofs (fe);

  const unsigned int level = tr.n_levels()-1;

  {
    deallog.push("ff");
    SparsityPattern bl(tr.n_cells(level-1), dof.n_dofs(level), (1<<dim) * fe.dofs_per_cell);;
    DoFTools::make_child_patches(bl, dof, level, false, false);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("tf");
    SparsityPattern bl(tr.n_cells(level-1), dof.n_dofs(level), (1<<dim) * fe.dofs_per_cell);;
    DoFTools::make_child_patches(bl, dof, level, true, false);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("tt");
    SparsityPattern bl(tr.n_cells(level-1), dof.n_dofs(level), (1<<dim) * fe.dofs_per_cell);;
    DoFTools::make_child_patches(bl, dof, level, true, true);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
}


int main()
{
  initlog();
  deallog.push("2D");
  test_global_refinement<2>(&test_block_list<2>);
  deallog.pop();
  deallog.push("3D");
  test_global_refinement<3>(&test_block_list<3>);
  deallog.pop();
}
