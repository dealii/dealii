// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check MGTransferMatrixFree with custom user constraints

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include "../tests.h"

template <int dim>
void
check()
{
  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr, 0, 1, false);

  typename Triangulation<dim>::active_cell_iterator cell = tr.begin_active();
  cell->face(0)->set_all_boundary_ids(1);

  tr.refine_global(1);

  FE_Q<dim> fe(1);

  DoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  mgdof.distribute_mg_dofs();

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(mgdof);

  const IndexSet relevant_dofs =
    DoFTools::extract_locally_relevant_level_dofs(mgdof, 0);
  AffineConstraints<double> user_constraints;
  user_constraints.reinit(mgdof.locally_owned_mg_dofs(0), relevant_dofs);

  typename DoFHandler<dim>::level_face_iterator face0 = mgdof.begin(0)->face(0);
  std::vector<types::global_dof_index>          face_dofs(fe.dofs_per_face);
  face0->get_mg_dof_indices(0, face_dofs);
  for (unsigned int i = 0; i < face_dofs.size(); ++i)
    {
      if (user_constraints.can_store_line(face_dofs[i]))
        {
          user_constraints.constrain_dof_to_zero(face_dofs[i]);
        }
    }
  user_constraints.close();
  mg_constrained_dofs.add_user_constraints(0, user_constraints);

  MGTransferMatrixFree<dim, double> transfer_mf(mg_constrained_dofs);
  transfer_mf.build(mgdof);

  deallog << "SRC Vector" << std::endl;
  LinearAlgebra::distributed::Vector<double> src_level_0(mgdof.n_dofs(0));
  src_level_0 = 1.0;
  for (unsigned int i = 0; i < mgdof.n_dofs(0); ++i)
    deallog << src_level_0(i) << ' ';
  deallog << std::endl << std::endl;

  {
    LinearAlgebra::distributed::Vector<double> dst_level_1(mgdof.n_dofs(1));
    transfer_mf.prolongate(1, dst_level_1, src_level_0);
    deallog << "DST Vector" << std::endl;
    for (unsigned int i = 0; i < mgdof.n_dofs(1); ++i)
      deallog << dst_level_1(i) << ' ';
    deallog << std::endl;
  }
}


int
main()
{
  initlog();

  check<2>();
  deallog << std::endl;
  check<3>();
  deallog << std::endl;
}
