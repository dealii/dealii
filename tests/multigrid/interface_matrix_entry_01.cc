// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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


// check functions MGTools::make_interface_sparsity_pattern()
// and MGConstrainedDoFs::is_interface_matrix_entry()

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tria, 0, 1);
  tria.refine_global(1);
  for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator
         cell = tria.begin_active();
       cell != tria.end();
       ++cell)
    for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
      {
        if (dim == 2)
          if (cell->vertex(v)[0] < 0.5 && cell->vertex(v)[1] < 0.5)
            {
              cell->set_refine_flag();
              break;
            }
        if (dim == 3)
          if (cell->vertex(v)[0] < 0.5 && cell->vertex(v)[1] < 0.5 &&
              cell->vertex(v)[2] < 0.5)
            {
              cell->set_refine_flag();
              break;
            }
      }
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>                 fe(2);
  DoFHandler<dim>           mg_dof_handler(tria);
  IndexSet                  locally_relevant_set;
  AffineConstraints<double> constraints;
  MGConstrainedDoFs         mg_constrained_dofs;

  mg_dof_handler.distribute_dofs(fe);
  mg_dof_handler.distribute_mg_dofs();

  DoFTools::extract_locally_relevant_dofs(mg_dof_handler, locally_relevant_set);

  constraints.reinit(locally_relevant_set);
  DoFTools::make_hanging_node_constraints(mg_dof_handler, constraints);

  std::set<types::boundary_id>                        dirichlet_boundary_ids;
  std::map<types::boundary_id, const Function<dim> *> dirichlet_boundary;
  Functions::ConstantFunction<dim> homogeneous_dirichlet_bc(0.0);
  dirichlet_boundary_ids.insert(0);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
  VectorTools::interpolate_boundary_values(mg_dof_handler,
                                           dirichlet_boundary,
                                           constraints);
  constraints.close();

  mg_constrained_dofs.clear();
  mg_constrained_dofs.initialize(mg_dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(mg_dof_handler,
                                                     dirichlet_boundary_ids);

  for (unsigned int level = 0; level < tria.n_levels(); ++level)
    {
      DynamicSparsityPattern dsp(mg_dof_handler.n_dofs(level),
                                 mg_dof_handler.n_dofs(level));
      MGTools::make_interface_sparsity_pattern(mg_dof_handler,
                                               mg_constrained_dofs,
                                               dsp,
                                               level);
      dsp.print(deallog.get_file_stream());
    }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
