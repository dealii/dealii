// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2015 by the deal.II authors
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


// check mg constrained dofs in parallel (different mesh than 01)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/function_map.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

template <int dim>
void setup_tria(parallel::distributed::Triangulation<dim> &triangulation)
{
  unsigned int n_subdiv = 1;
  GridGenerator::subdivided_hyper_cube (triangulation, n_subdiv, 0, 1);
  triangulation.refine_global(2);
  {
    for (typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active(); cell != triangulation.end(); ++cell)
      if (cell->is_locally_owned() &&
          cell->center().norm() < 0.55)
        cell->set_refine_flag();
    triangulation.execute_coarsening_and_refinement();
    for (typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active(); cell != triangulation.end(); ++cell)
      if (cell->is_locally_owned() &&
          cell->center().norm() > 0.3 && cell->center().norm() < 0.42)
        cell->set_refine_flag();
    triangulation.execute_coarsening_and_refinement();
    for (typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active(); cell != triangulation.end(); ++cell)
      if (cell->is_locally_owned() &&
          cell->center().norm() > 0.335 && cell->center().norm() < 0.39)
        cell->set_refine_flag();
    triangulation.execute_coarsening_and_refinement();
  }
}


template <typename DoFHandlerType>
void
extract_locally_active_level_dofs (const DoFHandlerType &dof_handler,
                                   const unsigned int    level,
                                   IndexSet             &dof_set)
{
  dof_set = IndexSet(dof_handler.n_dofs(level));

  // add all DoFs from our cells to the IndexSet

  // Note: For certain meshes (in particular in 3D and with many
  // processors), it is really necessary to cache intermediate data. After
  // trying several objects such as std::set, a vector that is always kept
  // sorted, and a vector that is initially unsorted and sorted once at the
  // end, the latter has been identified to provide the best performance.
  // Martin Kronbichler
  std::vector<types::global_dof_index> dof_indices;
  std::vector<types::global_dof_index> active_dofs;

  typename DoFHandlerType::cell_iterator cell = dof_handler.begin(level),
                                         endc = dof_handler.end(level);
  for (; cell!=endc; ++cell)
    {
      const types::subdomain_id id = cell->level_subdomain_id();
      if (id != dof_handler.get_triangulation().locally_owned_subdomain())
        continue;

      dof_indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_mg_dof_indices(dof_indices);
      for (unsigned int i=0; i<dof_indices.size(); ++i)
        if (!dof_set.is_element(dof_indices[i]))
          active_dofs.push_back(dof_indices[i]);
    }

  // sort, compress out duplicates, fill into index set
  std::sort(active_dofs.begin(), active_dofs.end());
  dof_set.add_indices(active_dofs.begin(), std::unique(active_dofs.begin(),
                                                       active_dofs.end()));

  dof_set.compress();
}




template <int dim>
void check_fe(FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
                                               Triangulation<dim>::none,
                                               parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  setup_tria(tr);

  ZeroFunction<dim> zero;
  typename FunctionMap<dim>::type fmap;
  fmap.insert(std::make_pair(0, &zero));

  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);
  dofh.distribute_mg_dofs(fe);

  MGConstrainedDoFs                    mg_constrained_dofs_ref;
  {
    // reorder
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_SELF,
                                                 Triangulation<dim>::none,
                                                 parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
    setup_tria(tr);

    DoFHandler<dim> dofhref(tr);
    dofhref.distribute_dofs(fe);
    dofhref.distribute_mg_dofs(fe);

    //std::map<std::string,std::vector<types::global_dof_index> > dofmap;
    std::map<std::string,std::vector<types::global_dof_index> > mgdofmap;

    for (typename DoFHandler<dim>::level_cell_iterator cell = dofhref.begin();
         cell != dofhref.end(); ++cell)
      {
        if (!cell->is_locally_owned_on_level())
          continue;

        std::vector<types::global_dof_index> &d = mgdofmap[cell->id().to_string()];
        d.resize(fe.dofs_per_cell);
        cell->get_mg_dof_indices(d);
      }

    for (typename DoFHandler<dim>::level_cell_iterator cell = dofh.begin();
         cell != dofh.end(); ++cell)
      {
        if (cell->level_subdomain_id()==numbers::artificial_subdomain_id)
          continue;

        std::vector<types::global_dof_index> &renumbered = mgdofmap[cell->id().to_string()];
        cell->set_mg_dof_indices(renumbered);
        cell->update_cell_dof_indices_cache();
      }

    typename FunctionMap<dim>::type      dirichlet_boundary;
    ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
    dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
    mg_constrained_dofs_ref.initialize(dofhref, dirichlet_boundary);
  }



  MGConstrainedDoFs                    mg_constrained_dofs;

  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
  mg_constrained_dofs.initialize(dofh, dirichlet_boundary);

  const unsigned int n_levels = tr.n_global_levels();
  for (unsigned int level = 0; level < n_levels; ++level)
    {
      deallog << "Level " << level << ":" << std::endl;

      IndexSet rei = mg_constrained_dofs.get_refinement_edge_indices (level);
      deallog << "get_refinement_edge_indices:" << std::endl;
      rei.print(deallog);

      IndexSet bi = mg_constrained_dofs.get_boundary_indices (level);
      deallog << "get_boundary_indices:" << std::endl;
      bi.print(deallog);

      IndexSet relevant;
      DoFTools::extract_locally_relevant_level_dofs (dofh,
                                                     level,
                                                     relevant);
      deallog << "relevant:" << std::endl;
      relevant.print(deallog);

      IndexSet active;
      extract_locally_active_level_dofs (dofh, level, active);
      deallog << "active:" << std::endl;
      active.print(deallog);

      // the indexsets should be the same when run in parallel (on the
      // active subset):
      deallog << (((rei & active) == (active & mg_constrained_dofs_ref.get_refinement_edge_indices(level)))
                  ?"ok ":"FAIL ")
              << (((bi & active) == (active & mg_constrained_dofs_ref.get_boundary_indices(level)))
                  ?"ok ":"FAIL ")
              << std::endl;

    }
}


template <int dim>
void check()
{
  FE_Q<dim> q1(1);
  FE_Q<dim> q2(2);
//  FE_DGQ<dim> dq1(1);

  FESystem<dim> s1(q1, 2, q2,1);

  check_fe(q1);
  //  check_fe(q2);
  //check_fe(s1);
}

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log (true);

  check<2> ();
  //check<3> ();
}
