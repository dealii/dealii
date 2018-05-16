// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// check that everything is ok when we have a triangulation that has
// fewer cells than there are processors
//
// this test is run with sufficiently many processors so that there
// are idle processors in 1d, 2d, and 3d
//
// this test is just like the one without hp_ but uses an
// hp::DoFHandler instead of a regular DoFHandler (but with only one
// element). the output should be, and is, the same


#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/base/utilities.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <numeric>


template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim>
  triangulation (MPI_COMM_WORLD, typename Triangulation<dim>::MeshSmoothing
                 (Triangulation<dim>::limit_level_difference_at_vertices), true,
                 parallel::shared::Triangulation<dim>::partition_custom_signal);
  triangulation.signals.post_refinement.connect ([&triangulation]()
  {
    // partition the triangulation by hand
    for (auto cell : triangulation.active_cell_iterators())
      cell->set_subdomain_id (cell->active_cell_index() %
                              Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD));
  }
                                                );

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (1);


  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));

  hp::DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);

  deallog << "n_dofs: " << dof_handler.n_dofs() << std::endl;
  deallog << "n_locally_owned_dofs: " << dof_handler.n_locally_owned_dofs() << std::endl;

  deallog << "n_locally_owned_dofs_per_processor: ";
  std::vector<types::global_dof_index> v = dof_handler.n_locally_owned_dofs_per_processor();
  unsigned int sum = 0;
  for (unsigned int i=0; i<v.size(); ++i)
    {
      deallog << v[i] << " ";
      sum += v[i];
    }
  deallog << " sum: " << sum << std::endl;
  deallog << " locally_owned_dofs: ";
  dof_handler.locally_owned_dofs().write(deallog.get_file_stream());
  deallog << std::endl;

  Assert (dof_handler.n_locally_owned_dofs() ==
          dof_handler.n_locally_owned_dofs_per_processor()[triangulation.locally_owned_subdomain()],
          ExcInternalError());
  Assert (dof_handler.n_locally_owned_dofs() == dof_handler.locally_owned_dofs().n_elements(),
          ExcInternalError());

  const unsigned int N = dof_handler.n_dofs();

  Assert (dof_handler.n_locally_owned_dofs() <= N,
          ExcInternalError());
  Assert (std::accumulate (dof_handler.n_locally_owned_dofs_per_processor().begin(),
                           dof_handler.n_locally_owned_dofs_per_processor().end(),
                           0U) == N,
          ExcInternalError());

  IndexSet all (N);
  for (unsigned int i=0;
       i<dof_handler.locally_owned_dofs_per_processor().size(); ++i)
    {
      IndexSet intersect = all & dof_handler.locally_owned_dofs_per_processor()[i];
      Assert(intersect.n_elements()==0, ExcInternalError());
      all.add_indices(dof_handler.locally_owned_dofs_per_processor()[i]);
    }

  Assert(all == complete_index_set(N), ExcInternalError());
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("1d");
  test<1>();
  deallog.pop();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
