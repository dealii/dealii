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


// Check clear_user_constraints()
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

  FE_Q<dim> fe = FE_Q<dim>(1);

  DoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  mgdof.distribute_mg_dofs();

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(mgdof);

  for (unsigned int level = 0; level < tr.n_levels(); ++level)
    {
      const IndexSet relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(mgdof, level);
      AffineConstraints<double> level_constraints;
      level_constraints.reinit(mgdof.locally_owned_mg_dofs(level),
                               relevant_dofs);

      typename DoFHandler<dim>::level_face_iterator face0 =
        mgdof.begin(level)->face(0);
      std::vector<types::global_dof_index> face_dofs(fe.dofs_per_face);
      face0->get_mg_dof_indices(level, face_dofs);
      for (unsigned int i = 0; i < face_dofs.size(); ++i)
        {
          if (level_constraints.can_store_line(face_dofs[i]))
            {
              level_constraints.constrain_dof_to_zero(face_dofs[i]);
              level_constraints.set_inhomogeneity(face_dofs[i], 5.0);
            }
        }
      level_constraints.close();
      mg_constrained_dofs.add_user_constraints(level, level_constraints);
      deallog << " level: " << level
              << " n_constraints: " << level_constraints.n_constraints()
              << std::endl;
    }

  mg_constrained_dofs.clear_user_constraints();
  for (unsigned int level = 0; level < tr.n_levels(); ++level)
    {
      const AffineConstraints<double> &level_constraints =
        mg_constrained_dofs.get_user_constraint_matrix(level);
      deallog << " level: " << level
              << " n_constraints: " << level_constraints.n_constraints()
              << std::endl;
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
