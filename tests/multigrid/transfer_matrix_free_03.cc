// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2016 by the deal.II authors
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


// Check MGTransferMatrixFree by comparison with MGTransferPrebuilt on a
// series of meshes with adaptive meshes for FE_Q (similar to
// transfer_matrix_free_02 but on a different mesh refining into a corner)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>


template <int dim, typename Number>
void check(const unsigned int fe_degree)
{
  deallog.threshold_double(std::max(5e2*(double)std::numeric_limits<Number>::epsilon(),
                                    1e-11));
  FE_Q<dim> fe(fe_degree);
  deallog << "FE: " << fe.get_name() << std::endl;

  // run a few different sizes...
  parallel::distributed::Triangulation<dim>
  tr(MPI_COMM_WORLD,
     Triangulation<dim>::limit_level_difference_at_vertices,
     parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube (tr, 3);
  for (unsigned int cycle=0; cycle<(dim == 2 ? 10 : 7); ++cycle)
    {
      // adaptive refinement into a circle
      for (typename Triangulation<dim>::active_cell_iterator cell=tr.begin_active(); cell != tr.end(); ++cell)
        if (cell->is_locally_owned() &&
            cell->vertex(0).norm() < 1e-10)
          cell->set_refine_flag();
      tr.execute_coarsening_and_refinement();

      deallog << "no. cells: " << tr.n_global_active_cells() << " on "
              << tr.n_global_levels() << " levels" << std::endl;

      DoFHandler<dim> mgdof(tr);
      mgdof.distribute_dofs(fe);
      mgdof.distribute_mg_dofs(fe);

      ConstraintMatrix hanging_node_constraints;
      IndexSet relevant_dofs;
      DoFTools::extract_locally_relevant_dofs(mgdof, relevant_dofs);
      hanging_node_constraints.reinit(relevant_dofs);
      DoFTools::make_hanging_node_constraints(mgdof, hanging_node_constraints);
      hanging_node_constraints.close();

      MGConstrainedDoFs mg_constrained_dofs;
      ZeroFunction<dim> zero_function;
      typename FunctionMap<dim>::type dirichlet_boundary;
      dirichlet_boundary[0] = &zero_function;
      mg_constrained_dofs.initialize(mgdof, dirichlet_boundary);

      // build reference
      MGTransferPrebuilt<parallel::distributed::Vector<double> >
      transfer_ref(hanging_node_constraints, mg_constrained_dofs);
      transfer_ref.build_matrices(mgdof);

      // build matrix-free transfer
      MGTransferMatrixFree<dim, Number> transfer(mg_constrained_dofs);
      transfer.build(mgdof);

      // check prolongation for all levels using random vector
      for (unsigned int level=1; level<mgdof.get_triangulation().n_global_levels(); ++level)
        {
          parallel::distributed::Vector<Number> v1, v2;
          parallel::distributed::Vector<double> v1_cpy, v2_cpy, v3;
          v1.reinit(mgdof.locally_owned_mg_dofs(level-1), MPI_COMM_WORLD);
          v2.reinit(mgdof.locally_owned_mg_dofs(level), MPI_COMM_WORLD);
          v3.reinit(mgdof.locally_owned_mg_dofs(level), MPI_COMM_WORLD);
          for (unsigned int i=0; i<v1.local_size(); ++i)
            v1.local_element(i) = (double)Testing::rand()/RAND_MAX;
          v1_cpy = v1;
          transfer.prolongate(level, v2, v1);
          transfer_ref.prolongate(level, v3, v1_cpy);
          v2_cpy = v2;
          v3 -= v2_cpy;
          deallog << "Diff prolongate   l" << level << ": " << v3.l2_norm() << std::endl;
        }

      // check restriction for all levels using random vector
      for (unsigned int level=1; level<mgdof.get_triangulation().n_global_levels(); ++level)
        {
          parallel::distributed::Vector<Number> v1, v2;
          parallel::distributed::Vector<double> v1_cpy, v2_cpy, v3;
          v1.reinit(mgdof.locally_owned_mg_dofs(level), MPI_COMM_WORLD);
          v2.reinit(mgdof.locally_owned_mg_dofs(level-1), MPI_COMM_WORLD);
          v3.reinit(mgdof.locally_owned_mg_dofs(level-1), MPI_COMM_WORLD);
          for (unsigned int i=0; i<v1.local_size(); ++i)
            v1.local_element(i) = (double)Testing::rand()/RAND_MAX;
          v1_cpy = v1;
          transfer.restrict_and_add(level, v2, v1);
          transfer_ref.restrict_and_add(level, v3, v1_cpy);
          v2_cpy = v2;
          v3 -= v2_cpy;
          deallog << "Diff restrict     l" << level << ": " << v3.l2_norm() << std::endl;

          v2 = 1.;
          v3 = 1.;
          transfer.restrict_and_add(level, v2, v1);
          transfer_ref.restrict_and_add(level, v3, v1_cpy);
          v2_cpy = v2;
          v3 -= v2_cpy;
          deallog << "Diff restrict add l" << level << ": " << v3.l2_norm() << std::endl;
        }
      deallog << std::endl;
    }
}


int main(int argc, char **argv)
{
  // no threading in this test...
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  mpi_initlog();

  check<2,double>(1);
  check<2,double>(3);
  check<3,double>(1);
  check<3,double>(3);
  check<2,float> (2);
  check<3,float> (2);
}
