// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check BlockSparsityPattern::write_svg

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_sparsity_pattern.h>

#include "../tests.h"

int
main()
{
  initlog();

  const int          dim = 2;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  const unsigned int degree = 1;
  FESystem<dim>      fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1);
  DoFHandler<dim>    dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::Cuthill_McKee(dof_handler);
  DoFRenumbering::component_wise(dof_handler);

  std::vector<unsigned int> block_component(dim + 1, 0);
  block_component[dim] = 1;
  const std::vector<types::global_dof_index> dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
  const unsigned int n_u = dofs_per_block[0];
  const unsigned int n_p = dofs_per_block[1];

  BlockSparsityPattern sparsity_pattern;
  {
    BlockDynamicSparsityPattern dsp(2, 2);

    dsp.block(0, 0).reinit(n_u, n_u);
    dsp.block(1, 0).reinit(n_p, n_u);
    dsp.block(0, 1).reinit(n_u, n_p);
    dsp.block(1, 1).reinit(n_p, n_p);

    dsp.collect_sizes();

    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);

    for (unsigned int c = 0; c < dim + 1; ++c)
      for (unsigned int d = 0; d < dim + 1; ++d)
        if (!((c == dim) && (d == dim)))
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;

    DoFTools::make_sparsity_pattern(dof_handler, coupling, dsp);

    sparsity_pattern.copy_from(dsp);
  }

  sparsity_pattern.print_svg(deallog.get_file_stream());
}
