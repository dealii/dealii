// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
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



// Test distributed solution transfer with fullydistributed triangulations.

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/solution_transfer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class InterpolationFunction : public Function<dim>
{
public:
  InterpolationFunction()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    return p.norm();
  }
};


using VectorType = LinearAlgebra::distributed::Vector<double>;


template <int dim>
VectorType
transfer(const std::string &mode)
{
  AssertThrow(mode == "save" || mode == "load", ExcInternalError());

  const std::string filename =
    "solution_transfer_" + std::to_string(dim) + "d_out";

  parallel::fullydistributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  {
    parallel::distributed::Triangulation<dim> tria_base(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(tria_base);
    tria_base.refine_global(3);

    const auto description = TriangulationDescription::Utilities::
      create_description_from_triangulation(tria_base, MPI_COMM_WORLD);

    triangulation.create_triangulation(description);
  }

  if (mode == "load")
    {
      triangulation.clear();
      triangulation.load(filename);
    }

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(FE_Q<dim>(2));

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
  std::shared_ptr<Utilities::MPI::Partitioner> partitioner =
    std::make_shared<Utilities::MPI::Partitioner>(
      dof_handler.locally_owned_dofs(), locally_relevant_dofs, MPI_COMM_WORLD);

  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer(
    dof_handler);
  VectorType vector(partitioner);
  if (mode == "save")
    {
      VectorTools::interpolate(dof_handler,
                               InterpolationFunction<dim>(),
                               vector);
      vector.update_ghost_values();

      solution_transfer.prepare_for_serialization(vector);

      triangulation.save(filename);
    }
  else if (mode == "load")
    {
      solution_transfer.deserialize(vector);

      vector.update_ghost_values();
    }

  return vector;
}


template <int dim>
void
test()
{
  VectorType vector        = transfer<dim>("save");
  VectorType vector_loaded = transfer<dim>("load");

  // Verify that both vectors are the same.
  vector.add(-1, vector_loaded);
  deallog << (vector.linfty_norm() < 1e-16 ? "PASSED" : "FAILED") << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
