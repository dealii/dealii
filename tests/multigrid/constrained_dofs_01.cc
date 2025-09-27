// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check mg constrained dofs in parallel

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
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>

#include <algorithm>

#include "../tests.h"



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
      if (cell->id().to_string() == "0_2:11")
        cell->set_refine_flag();
    }
  tr.execute_coarsening_and_refinement();
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

  Functions::ZeroFunction<dim>                        zero;
  std::map<types::boundary_id, const Function<dim> *> fmap;
  fmap.insert(std::make_pair(0, &zero));

  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);
  dofh.distribute_mg_dofs();

  MGConstrainedDoFs mg_constrained_dofs_ref;
  {
    // reorder
    parallel::distributed::Triangulation<dim> tr(
      MPI_COMM_SELF,
      Triangulation<dim>::limit_level_difference_at_vertices,
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
    setup_tria(tr);

    DoFHandler<dim> dofhref(tr);
    dofhref.distribute_dofs(fe);
    dofhref.distribute_mg_dofs();

    // std::map<std::string,std::vector<types::global_dof_index> > dofmap;
    std::map<std::string, std::vector<types::global_dof_index>> mgdofmap;

    for (typename DoFHandler<dim>::level_cell_iterator cell = dofhref.begin();
         cell != dofhref.end();
         ++cell)
      {
        if (!cell->is_locally_owned_on_level())
          continue;

        std::vector<types::global_dof_index> &d =
          mgdofmap[cell->id().to_string()];
        d.resize(fe.dofs_per_cell);
        cell->get_mg_dof_indices(d);
      }

    for (typename DoFHandler<dim>::level_cell_iterator cell = dofh.begin();
         cell != dofh.end();
         ++cell)
      {
        if (cell->level_subdomain_id() == numbers::artificial_subdomain_id)
          continue;

        std::vector<types::global_dof_index> &renumbered =
          mgdofmap[cell->id().to_string()];
        cell->set_mg_dof_indices(renumbered);
      }

    mg_constrained_dofs_ref.initialize(dofhref);
    mg_constrained_dofs_ref.make_zero_boundary_constraints(dofhref, {0});
  }



  MGConstrainedDoFs mg_constrained_dofs;

  mg_constrained_dofs.initialize(dofh);
  mg_constrained_dofs.make_zero_boundary_constraints(dofh, {0});

  const unsigned int n_levels = tr.n_global_levels();
  for (unsigned int level = 0; level < n_levels; ++level)
    {
      deallog << "Level " << level << ':' << std::endl;

      IndexSet rei = mg_constrained_dofs.get_refinement_edge_indices(level);
      deallog << "get_refinement_edge_indices:" << std::endl;
      rei.print(deallog);

      IndexSet bi = mg_constrained_dofs.get_boundary_indices(level);
      deallog << "get_boundary_indices:" << std::endl;
      bi.print(deallog);

      const IndexSet relevant =
        DoFTools::extract_locally_relevant_level_dofs(dofh, level);
      deallog << "relevant:" << std::endl;
      relevant.print(deallog);

      // the indexsets should be the same when run in parallel (on the
      // relevant subset):
      deallog
        << ((rei ==
             (relevant &
              mg_constrained_dofs_ref.get_refinement_edge_indices(level))) ?
              "ok " :
              "FAIL ")
        << ((bi ==
             (relevant & mg_constrained_dofs_ref.get_boundary_indices(level))) ?
              "ok " :
              "FAIL ")
        << std::endl;
    }

  {
    deallog << "Test add_boundary_indices: " << std::endl;
    // used to test add_boundary_indices()
    MGConstrainedDoFs mg_constrained_dofs_1;
    mg_constrained_dofs_1.initialize(dofh);

    MGLevelObject<AffineConstraints<double>> mg_boundary_constraints;
    mg_boundary_constraints.resize(0, n_levels - 1);

    std::vector<IndexSet> boundary_indices;
    boundary_indices.resize(n_levels);
    MGTools::make_boundary_list(dofh, {0}, boundary_indices);

    for (unsigned int level = 0; level < n_levels; ++level)
      {
        deallog << "Level " << level << ':' << std::endl;

        const IndexSet relevant =
          DoFTools::extract_locally_relevant_level_dofs(dofh, level);
        mg_boundary_constraints[level].reinit(dofh.locally_owned_mg_dofs(level),
                                              relevant);
        mg_boundary_constraints[level].add_lines(boundary_indices[level]);

        mg_constrained_dofs_1.add_boundary_indices(dofh,
                                                   level,
                                                   boundary_indices[level]);
        const auto &bi = mg_constrained_dofs_1.get_boundary_indices(level);

        deallog << "get_boundary_indices test:" << std::endl;
        bi.print(deallog);

        deallog << "relevant:" << std::endl;
        relevant.print(deallog);

        deallog << ((bi ==
                     (relevant &
                      mg_constrained_dofs_ref.get_boundary_indices(level))) ?
                      "ok " :
                      "FAIL test")
                << std::endl;
      }

    // extract boundary indices from constraint matrices
    // this is probably how we would use the function add_boundary_indices()
    mg_constrained_dofs_1.clear();
    mg_constrained_dofs_1.initialize(dofh);
    for (unsigned int level = 0; level < n_levels; ++level)
      {
        deallog << "Level " << level << ':' << std::endl;

        IndexSet level_boundary_indices(dofh.n_dofs(level));
        for (auto line : mg_boundary_constraints[level].get_lines())
          {
            if (line.entries.size() == 0)
              {
                level_boundary_indices.add_index(line.index);
              }
          }
        mg_constrained_dofs_1.add_boundary_indices(dofh,
                                                   level,
                                                   level_boundary_indices);
        const auto &bi = mg_constrained_dofs_1.get_boundary_indices(level);

        const IndexSet relevant =
          DoFTools::extract_locally_relevant_level_dofs(dofh, level);

        deallog << ((bi ==
                     (relevant &
                      mg_constrained_dofs_ref.get_boundary_indices(level))) ?
                      "ok " :
                      "FAIL test")
                << std::endl;
      }
  }
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
