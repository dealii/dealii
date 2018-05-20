// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// tests DataOutFaces in parallel for a case where one processor has
// no faces to deal with

// This test case shows that the DataOutFaces class is able to
// support parallel triangulations.  The main problem that I was
// experiencing was the mesh I was testing on was too coarse for
// larger number of processors. This test case shows that as
// well. For 4 processors the code produces output without error
// for both the 12 repitions and the 2 repititions. For 6 and 12
// processors only the 12 repetition case produces the proper
// output. Fortunately it does show as long as the mesh is
// adequately refined DataOutFaces produces the output for each
// subdomain.

#include "../tests.h"
#include <deal.II/lac/generic_linear_algebra.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/numerics/data_out_faces.h>

int
main(int argc, char* argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  {
    const unsigned int dim = 2;

    MPI_Comm                                  mpi_communicator(MPI_COMM_WORLD);
    parallel::distributed::Triangulation<dim> triangulation(
      mpi_communicator,
      typename Triangulation<dim>::MeshSmoothing(
        Triangulation<dim>::smoothing_on_refinement
        | Triangulation<dim>::smoothing_on_coarsening));
    FESystem<dim>   fe(FE_FaceQ<dim>(2), dim);
    DoFHandler<dim> dof_handler(triangulation);

    triangulation.clear();
    GridGenerator::subdivided_hyper_cube(triangulation, 1, 0, 1);
    dof_handler.distribute_dofs(fe);

    LinearAlgebraTrilinos::MPI::Vector locally_relevant_sol;
    IndexSet                           locally_owned_dofs;
    IndexSet                           locally_relevant_dofs;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    locally_relevant_sol.reinit(
      locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    DataOutFaces<dim> data_out_face(false);

    data_out_face.add_data_vector(dof_handler, locally_relevant_sol, "data");

    data_out_face.build_patches();
    data_out_face.write_gnuplot(deallog.get_file_stream());
  }
}
