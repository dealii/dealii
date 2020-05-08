// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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



// Check whether chains of constraints traverse over dofs on faces,
// and if this is an issue in parallel applications.
//
// Domain in 2D/3D:
//           +---------+
//           |         |
//           | FE_Q(1) |
//           |         |
// +---------+---------+---------+
// |         |         |         |
// | FE_Q(1) | FE_Q(2) | FE_Q(1) |
// |         |         |         |
// +---------+---------+---------+


#include <deal.II/base/geometry_info.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree_center,
     const unsigned int degree_other,
     const bool         print_constraints = false)
{
  Assert(dim > 1, ExcNotImplemented());

  // ------ setup ------
  // prepare Triangulation as a hyper_cross
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  std::vector<unsigned int> sizes;
  sizes.insert(sizes.end(), {1, 1}); // extension in in x-direction
  sizes.insert(sizes.end(), {0, 1}); // extension in in y-direction
  if (dim > 2)
    sizes.insert(sizes.end(), {0, 0}); // extension in in z-direction

  GridGenerator::hyper_cross(tria, sizes);

  // prepare FECollection with arbitrary number of entries
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FESystem<dim>(FE_Q<dim>(degree_other), dim));
  fe_collection.push_back(FESystem<dim>(FE_Q<dim>(degree_center), dim));

  // prepare DoFHandler
  hp::DoFHandler<dim> dh(tria);

  for (const auto &cell : dh.active_cell_iterators())
    if (cell->is_locally_owned() && cell->id().to_string() == "1_0:")
      {
        // set different fe on center cell
        cell->set_active_fe_index(1);

#ifdef DEBUG
        // verify that our scenario is initialized correctly
        // by checking the number of neighbors of the center cell
        unsigned int n_neighbors = 0;
        for (const unsigned int i : GeometryInfo<dim>::face_indices())
          if (static_cast<unsigned int>(cell->neighbor_index(i)) !=
              numbers::invalid_unsigned_int)
            ++n_neighbors;
        Assert(n_neighbors == 3, ExcInternalError());
#endif
      }

  dh.distribute_dofs(fe_collection);

  // ---- constrain ----
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dh, locally_relevant_dofs);

  AffineConstraints<double> constraints;
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(dh, constraints);

  constraints.close();

  // ------ verify -----
  IndexSet locally_active_dofs;
  DoFTools::extract_locally_active_dofs(dh, locally_active_dofs);

  if (print_constraints)
    {
      deallog << "constraints:" << std::endl;
      constraints.print(deallog.get_file_stream());
    }
  deallog << "consistent? "
          << constraints.is_consistent_in_parallel(
               dh.locally_owned_dofs_per_processor(),
               locally_active_dofs,
               MPI_COMM_WORLD,
               true)
          << std::endl;

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>(2, 1);
  deallog.pop();

  deallog.push("3d");
  test<3>(2, 1);
  deallog.pop();
}
