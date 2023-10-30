// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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


/**
 * Compare no-normal-flux constrains in the case of MGTransferMatrixFree
 * and MGTransferMF.
 */

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria(
    Triangulation<dim>::MeshSmoothing::limit_level_difference_at_vertices);
  GridGenerator::hyper_ball(tria);
  tria.refine_global(1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(FE_Q<dim>(1), dim));
  dof_handler.distribute_mg_dofs();

  MGConstrainedDoFs mg_constrained_dofs;

  mg_constrained_dofs.initialize(dof_handler);

  if (false)
    {
      mg_constrained_dofs.make_no_normal_flux_constraints(dof_handler, 0, 0);
    }
  else
    {
      for (unsigned int level = 0; level < 2; ++level)
        {
          AffineConstraints<double> user_level_constraints;
          user_level_constraints.reinit(
            dof_handler.locally_owned_mg_dofs(level));

          const IndexSet &refinement_edge_indices =
            mg_constrained_dofs.get_refinement_edge_indices(level);

          const MappingQ1<dim> mapping;

          VectorTools::compute_no_normal_flux_constraints_on_level(
            dof_handler,
            0,
            {0},
            user_level_constraints,
            mapping,
            refinement_edge_indices,
            level);
          user_level_constraints.close();

          mg_constrained_dofs.add_user_constraints(level,
                                                   user_level_constraints);

          user_level_constraints.print(deallog.get_file_stream());
          deallog << std::endl;
        }
    }

  LinearAlgebra::distributed::Vector<double> vec_coarse(dof_handler.n_dofs(0));
  LinearAlgebra::distributed::Vector<double> vec_fine(dof_handler.n_dofs(1));

  // version 1: MGTransferMatrixFree
  MGTransferMatrixFree<dim, double> transfer_0(mg_constrained_dofs);
  transfer_0.build(dof_handler);

  vec_coarse = 1.0;
  vec_fine   = 0.0;
  transfer_0.prolongate(1, vec_fine, vec_coarse);
  vec_fine.print(deallog.get_file_stream());
  deallog << std::endl;

  vec_coarse = 0.0;
  vec_fine   = 1.0;
  transfer_0.restrict_and_add(1, vec_coarse, vec_fine);
  mg_constrained_dofs.get_user_constraint_matrix(0).set_zero(vec_coarse);
  vec_coarse.print(deallog.get_file_stream());
  deallog << std::endl;

  // version 2: MGTransferMF
  MGTransferMF<dim, double> transfer_1(mg_constrained_dofs);
  transfer_1.build(dof_handler);

  vec_coarse = 1.0;
  vec_fine   = 0.0;
  transfer_1.prolongate(1, vec_fine, vec_coarse);
  vec_fine.print(deallog.get_file_stream());
  deallog << std::endl;

  vec_coarse = 0.0;
  vec_fine   = 1.0;
  transfer_1.restrict_and_add(1, vec_coarse, vec_fine);
  mg_constrained_dofs.get_user_constraint_matrix(0).set_zero(vec_coarse);
  vec_coarse.print(deallog.get_file_stream());
  deallog << std::endl;
}

int
main()
{
  initlog();

  test<2>();
}
