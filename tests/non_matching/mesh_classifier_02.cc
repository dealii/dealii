// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check that the MeshClassifier throws an exception in debug mode if the
// ghost entries of a LinearAlgebra::distributed::Vector do not cover all
// locally relevant DoFs, as is the case for vectors initialized by MatrixFree,
// and that it works if all locally relevant DoFs are present.

#include <deal.II/base/function_signed_distance.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/non_matching/mesh_classifier.h>

#include <deal.II/numerics/vector_tools_interpolate.h>

#include "../tests.h"


template <int dim>
void
test(const bool use_locally_relevant_dofs)
{
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(3);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  // The locally active DoFs correspond to the ghost entries of a vector
  // initialized by MatrixFree.
  const IndexSet ghost_dofs =
    use_locally_relevant_dofs ?
      DoFTools::extract_locally_relevant_dofs(dof_handler) :
      DoFTools::extract_locally_active_dofs(dof_handler);

  LinearAlgebra::distributed::Vector<double> level_set(
    dof_handler.locally_owned_dofs(), ghost_dofs, MPI_COMM_WORLD);

  const Functions::SignedDistance::Sphere<dim> sphere(Point<dim>(), 0.5);
  VectorTools::interpolate(dof_handler, sphere, level_set);
  level_set.update_ghost_values();

  try
    {
      NonMatching::MeshClassifier<dim> classifier(dof_handler, level_set);
      classifier.reclassify();
      deallog << "Classification succeeded" << std::endl;
    }
  catch (const ExceptionBase &)
    {
      deallog << "Exception caught" << std::endl;
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  MPILogInitAll                    log;
  deal_II_exceptions::disable_abort_on_exception();

  deallog << "locally active DoFs as ghosts" << std::endl;
  test<2>(false);

  deallog << "locally relevant DoFs as ghosts" << std::endl;
  test<2>(true);
}
