// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2016 by the deal.II authors
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


// Deadlock reported by Kronbichler (github
// https://github.com/dealii/dealii/issues/2051) with 3 processes
// in MgTransferPrebuilt


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/grid/grid_out.h>


template <int dim>
void check()
{
  FE_Q<dim> fe(1);

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

      unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

      if (0)
        {
          std::ofstream grid_output (("out"+Utilities::to_string(myid)+".svg").c_str());
          GridOut grid_out;
          GridOutFlags::Svg flags;
          flags.label_level_subdomain_id = true;
          flags.coloring = GridOutFlags::Svg::level_subdomain_id;
          flags.convert_level_number_to_height = true;
          grid_out.set_flags(flags);

          grid_out.write_svg (tr, grid_output);
        }

      MGTransferPrebuilt<LinearAlgebra::distributed::Vector<double> >
      transfer_ref(hanging_node_constraints, mg_constrained_dofs);
      transfer_ref.build_matrices(mgdof);
    }
}

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  mpi_initlog();

  check<2>();
}
