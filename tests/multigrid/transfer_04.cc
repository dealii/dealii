// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


// check mg transfer in parallel, especially communication of copy_indices

#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>

#include <deal.II/numerics/data_out.h>

#include <algorithm>

#include "../tests.h"

using namespace std;

template <int dim>
void
setup_tria(parallel::distributed::Triangulation<dim> &tr)
{
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator
         cell = tr.begin_active();
       cell != tr.end();
       ++cell)
    {
      if (cell->id().to_string() == "0_2:03" ||
          cell->id().to_string() == "0_2:00" ||
          cell->id().to_string() == "0_2:01" ||
          cell->id().to_string() == "0_2:12")
        cell->set_refine_flag();
    }
  tr.execute_coarsening_and_refinement();
  for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator
         cell = tr.begin_active();
       cell != tr.end();
       ++cell)
    {
      if (cell->id().to_string() == "0_3:032" ||
          cell->id().to_string() == "0_3:000")
        cell->set_refine_flag();
    }
  tr.execute_coarsening_and_refinement();


  for (typename parallel::distributed::Triangulation<dim>::cell_iterator cell =
         tr.begin();
       cell != tr.end();
       ++cell)
    {
      deallog << "cell=" << cell->id()
              << " level_subdomain_id=" << cell->level_subdomain_id()
              << std::endl;
    }
}



template <int dim>
void
check_fe(FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  setup_tria(tr);

  if (false)
    {
      DataOut<dim>  data_out;
      Vector<float> subdomain(tr.n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = tr.locally_owned_subdomain();
      data_out.attach_triangulation(tr);
      data_out.add_data_vector(subdomain, "subdomain");
      data_out.build_patches(0);
      const std::string filename =
        ("solution." +
         Utilities::int_to_string(tr.locally_owned_subdomain(), 4) + ".vtu");
      std::ofstream output(filename.c_str());
      data_out.write_vtu(output);
    }

  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);
  dofh.distribute_mg_dofs();
  typedef TrilinosWrappers::MPI::Vector vector_t;

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dofh);

  AffineConstraints<double> hanging_node_constraints;
  IndexSet                  locally_relevant_set;
  DoFTools::extract_locally_relevant_dofs(dofh, locally_relevant_set);
  hanging_node_constraints.reinit(locally_relevant_set);
  DoFTools::make_hanging_node_constraints(dofh, hanging_node_constraints);
  hanging_node_constraints.close();

  MGTransferPrebuilt<vector_t> transfer(mg_constrained_dofs);
  transfer.build(dofh);
  // transfer.print_indices(deallog.get_file_stream());

  MGLevelObject<vector_t> u(0, tr.n_global_levels() - 1);
  for (unsigned int level = u.min_level(); level <= u.max_level(); ++level)
    {
      u[level].reinit(dofh.locally_owned_mg_dofs(level), MPI_COMM_WORLD);
      for (unsigned int i = 0;
           i < dofh.locally_owned_mg_dofs(level).n_elements();
           ++i)
        {
          unsigned int index =
            dofh.locally_owned_mg_dofs(level).nth_index_in_set(i);
          u[level][index] = 1.0; // 1000+level*100+index;
        }
      u[level].compress(VectorOperation::insert);
    }

  vector_t v;
  v.reinit(dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  v = 0.;
  transfer.copy_from_mg(dofh, v, u);
  hanging_node_constraints.distribute(v);

  {
    for (unsigned int i = 0; i < dofh.locally_owned_dofs().n_elements(); ++i)
      {
        unsigned int index = dofh.locally_owned_dofs().nth_index_in_set(i);
        if (std::abs(v[index] - 1.0) > 1e-5)
          deallog << "ERROR: index=" << index << " is equal to " << v[index]
                  << std::endl;
      }
  }
  deallog << "ok" << std::endl;
}


template <int dim>
void
check()
{
  FE_Q<dim> q1(1);
  FE_Q<dim> q2(2);
  //  FE_DGQ<dim> dq1(1);

  FESystem<dim> s1(q1, 2, q2, 1);

  check_fe(q1);
  //  check_fe(q2);
  // check_fe(s1);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  check<2>();
  // check<3> ();
}
