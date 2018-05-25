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



// have a 2x1 coarse mesh (or 2x1x1) and verify DoF indices in the hp
// case with an FECollection that contains multiple copies of the same
// FE_Q(2) element. in the sequential case, the hp code will unify DoF
// indices on boundaries between locally owned subdomains, but in
// parallel, we don't do that


#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>

#include <numeric>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices);

  std::vector<unsigned int> reps(dim, 1U);
  reps[0] = 2;
  Point<dim> top_right;
  for (unsigned int d = 0; d < dim; ++d)
    top_right[d] = (d == 0 ? 2 : 1);
  GridGenerator::subdivided_hyper_rectangle(
    triangulation, reps, Point<dim>(), top_right);
  Assert(triangulation.n_global_active_cells() == 2, ExcInternalError());
  Assert(triangulation.n_active_cells() == 2, ExcInternalError());

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(2));
  fe.push_back(FE_Q<dim>(2));

  hp::DoFHandler<dim> dof_handler(triangulation);
  if (dof_handler.begin_active()->is_locally_owned())
    dof_handler.begin_active()->set_active_fe_index(0);
  if ((++dof_handler.begin_active())->is_locally_owned())
    (++dof_handler.begin_active())->set_active_fe_index(1);
  dof_handler.distribute_dofs(fe);

  deallog << "Processor: " << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
          << std::endl;
  for (auto cell : dof_handler.active_cell_iterators())
    {
      deallog << "  Cell: " << cell << std::endl;

      std::vector<types::global_dof_index> dof_indices(
        cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(dof_indices);
      deallog << "    ";
      for (auto i : dof_indices)
        deallog << i << ' ';
      deallog << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
