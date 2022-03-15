// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

// Test DataOutResample to create a slice.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/data_out_resample.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
class AnalyticalFunction : public Function<dim>
{
public:
  AnalyticalFunction()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    (void)component;

    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
  }
};


template <int dim, int spacedim>
std::shared_ptr<const Utilities::MPI::Partitioner>
create_partitioner(const DoFHandler<dim, spacedim> &dof_handler)
{
  IndexSet locally_relevant_dofs;

  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  return std::make_shared<const Utilities::MPI::Partitioner>(
    dof_handler.locally_owned_dofs(),
    locally_relevant_dofs,
    dof_handler.get_communicator());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  const int dim       = 3;
  const int patch_dim = 2;
  const int spacedim  = 3;

  const unsigned int n_refinements_1 = 3;
  const unsigned int n_refinements_2 = 3;
  const MPI_Comm     comm            = MPI_COMM_WORLD;

  parallel::distributed::Triangulation<patch_dim, spacedim> tria_slice(comm);
  GridGenerator::hyper_cube(tria_slice, -1.0, +1.0);
  tria_slice.refine_global(n_refinements_2);

  MappingQ1<patch_dim, spacedim> mapping_slice;

  parallel::distributed::Triangulation<dim, spacedim> tria_backround(comm);
  GridGenerator::hyper_cube(tria_backround, -1.0, +1.0);
  tria_backround.refine_global(n_refinements_1);

  DoFHandler<dim, spacedim> dof_handler(tria_backround);
  dof_handler.distribute_dofs(FE_Q<dim, spacedim>{1});

  MappingQ1<dim, spacedim> mapping;

  LinearAlgebra::distributed::Vector<double> vector(
    create_partitioner(dof_handler));

  VectorTools::interpolate(mapping,
                           dof_handler,
                           AnalyticalFunction<dim>(),
                           vector);

  vector.update_ghost_values();

  DataOutResample<dim, patch_dim, spacedim> data_out(tria_slice, mapping_slice);
  data_out.add_data_vector(dof_handler, vector, "solution_0");
  data_out.add_data_vector(dof_handler, vector, "solution_1");
  data_out.update_mapping(mapping);
  data_out.build_patches();

#if 1
  data_out.write_vtk(deallog.get_file_stream());
#else
  data_out.write_vtu_with_pvtu_record("./", "data_out_01", 0, comm, 1, 1);
#endif
}
