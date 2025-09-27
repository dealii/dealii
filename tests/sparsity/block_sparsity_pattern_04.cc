// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Verify that AffineConstraints::add_entries_local_to_global() works correctly.
// In particular, there used to be (written before
// BlockSparsityPattern::add_entries() existed) a block-specific version of that
// function. This test verifies that the non-block version, combined with the
// block version of add_entries(), produces the same output.

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include <deal.II/numerics/vector_tools_boundary.h>

#include <vector>

#include "../tests.h"

int
main()
{
  initlog();

  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FESystem<2>        fe(FE_Q<2>(2), 2, FE_DGQ<2>(1), 1, FE_DGQ<2>(0), 2);
  const unsigned int n_components = fe.n_components();
  const unsigned int n_blocks     = fe.n_blocks();
  DoFHandler<2>      dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  VectorTools::interpolate_boundary_values(
    MappingCartesian<2>(),
    dof_handler,
    0,
    Functions::ConstantFunction<2>(42.0, n_components),
    constraints,
    ComponentMask({true, true, false, false, false}));
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  Table<2, DoFTools::Coupling> coupling(n_components, n_components);
  for (unsigned int c = 0; c < n_components; ++c)
    for (unsigned int d = 0; d < n_components; ++d)
      // This is completely arbitrary for this test so just flip back and forth
      coupling[c][d] = (c + d) % 2 == 0 ? DoFTools::always : DoFTools::none;

  const std::vector<types::global_dof_index> dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler);

  BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
  DoFTools::make_sparsity_pattern(
    dof_handler, coupling, dsp, constraints, false);

  BlockSparsityPattern sp(n_blocks, n_blocks);
  sp.copy_from(dsp);
  sp.print(deallog.get_file_stream());
}
