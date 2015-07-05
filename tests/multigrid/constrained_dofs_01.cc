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


// check mg constrained dofs in parallel

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

  std::string id_to_string(const CellId &id)
  {
    std::ostringstream ss;
    ss << id;
    return ss.str();
  }

template <int dim>
void setup_tria(parallel::distributed::Triangulation<dim> &tr)
{
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  
  for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = tr.begin_active();
       cell != tr.end(); ++cell)
    {
      if (id_to_string(cell->id()) == "0_2:11")
	cell->set_refine_flag();
    }
  tr.execute_coarsening_and_refinement();
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
  { // reorder
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

      std::vector<types::global_dof_index> &d = mgdofmap[id_to_string(cell->id())];
      d.resize(fe.dofs_per_cell);
      cell->get_mg_dof_indices(d);
    }

  for (typename DoFHandler<dim>::level_cell_iterator cell = dofh.begin();
       cell != dofh.end(); ++cell)
    {
      if (cell->level_subdomain_id()==numbers::artificial_subdomain_id)
	continue;
      
      std::vector<types::global_dof_index> &renumbered = mgdofmap[id_to_string(cell->id())];
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
      DoFTools::extract_locally_relevant_mg_dofs (dofh,
                relevant, level);
      deallog << "relevant:" << std::endl;
      relevant.print(deallog);

      // the indexsets should be the same when run in parallel (on the
      // relevant subset):
      deallog << ((rei == (relevant & mg_constrained_dofs_ref.get_refinement_edge_indices(level)))
		  ?"ok ":"FAIL ")
	      << ((bi == (relevant & mg_constrained_dofs_ref.get_boundary_indices(level)))
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
  MPILogInitAll log;

  check<2> ();
  //check<3> ();
}
