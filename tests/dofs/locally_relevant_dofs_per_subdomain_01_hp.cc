// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2018 by the deal.II authors
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


// Test DoFTools::locally_relevant_dofs_per_subdomain
// using a refined shared triangulation
//
// like the test without _hp, but for a hp::DoFHandler

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
void
write_mesh(const parallel::shared::Triangulation<dim> &tria,
           const char *                                filename_)
{
  DataOut<dim> data_out;
  data_out.attach_triangulation(tria);
  Vector<float> subdomain(tria.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = tria.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");

  data_out.build_patches();
  const std::string filename =
    (filename_ + Utilities::int_to_string(tria.locally_owned_subdomain(), 4));
  {
    std::ofstream output((filename + ".vtu").c_str());
    data_out.write_vtu(output);
  }
}



template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_zorder);

  hp::FECollection<dim> fe;
  fe.push_back(FESystem<dim>(FE_Q<dim>(3), 2, FE_DGQ<dim>(1), 1));

  DoFHandler<dim> dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  dof_handler.distribute_dofs(fe);

  // write_mesh(tr, "mesh");

  const std::vector<IndexSet> locally_relevant_dofs_per_subdomain =
    DoFTools::locally_relevant_dofs_per_subdomain(dof_handler);

  deallog << "locally_relevant_dofs on subdomain "
          << triangulation.locally_owned_subdomain() << ": ";
  locally_relevant_dofs_per_subdomain[triangulation.locally_owned_subdomain()]
    .print(deallog);
  deallog << "\n" << std::endl;
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
}
