// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
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
  dof.distribute_mg_dofs();

  const unsigned int level = tr.n_levels() - 1;

  {
    deallog.push("t");
    SparsityPattern bl;
    DoFTools::make_single_patch(bl, dof, level, true);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("f");
    SparsityPattern bl;
    DoFTools::make_single_patch(bl, dof, level, false);
    bl.compress();
    print_patches(bl);
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
