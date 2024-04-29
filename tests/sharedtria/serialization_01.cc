// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test serialization with shared triangulations.

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "../grid/tests.h"


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

template <int dim, typename TriangulationType>
void
test(TriangulationType &triangulation, TriangulationType &triangulation_clean)
{
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(FE_Q<dim>(2));

  using VectorType = Vector<double>;

  VectorType vector(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, InterpolationFunction<dim>(), vector);

  VectorType vector_loaded(dof_handler.n_dofs());

  print_statistics(triangulation, false);

  std::stringstream stream_0, stream_1;
  {
    boost::archive::text_oarchive oa_0(stream_0);
    boost::archive::text_oarchive oa_1(stream_1);
    oa_0 << triangulation;
    oa_0 << vector;
    oa_1 << triangulation;
    oa_1 << vector;
  }

  triangulation.clear();

  {
    // load into existing tringulation
    {
      boost::archive::text_iarchive ia(stream_0);
      ia >> triangulation;
      ia >> vector_loaded;
    }
    print_statistics(triangulation, false);

    // Verify that error is 0.
    VectorType error(vector);
    error.add(-1, vector_loaded);

    deallog << (error.linfty_norm() < 1e-16 ? "PASSED" : "FAILED") << std::endl;
  }

  {
    // load into new tringulation
    {
      boost::archive::text_iarchive ia(stream_1);
      ia >> triangulation_clean;
      ia >> vector_loaded;
    }
    print_statistics(triangulation_clean, false);

    // Verify that error is 0.
    VectorType error(vector);
    error.add(-1, vector_loaded);

    deallog << (error.linfty_norm() < 1e-16 ? "PASSED" : "FAILED") << std::endl;
  }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.push("2d");
  {
    constexpr int dim = 2;

    parallel::shared::Triangulation<dim> triangulation(
      MPI_COMM_WORLD,
      ::Triangulation<dim>::none,
      false,
      parallel::shared::Triangulation<dim>::partition_zorder);
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(3);


    parallel::shared::Triangulation<dim> triangulation_other(
      MPI_COMM_WORLD,
      ::Triangulation<dim>::none,
      false,
      parallel::shared::Triangulation<dim>::partition_zorder);

    test<dim>(triangulation, triangulation_other);
  }
  deallog.pop();

  deallog.push("3d");
  {
    constexpr int dim = 3;

    parallel::shared::Triangulation<dim> triangulation(
      MPI_COMM_WORLD,
      ::Triangulation<dim>::none,
      false,
      parallel::shared::Triangulation<dim>::partition_zorder);
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(3);

    parallel::shared::Triangulation<dim> triangulation_other(
      MPI_COMM_WORLD,
      ::Triangulation<dim>::none,
      false,
      parallel::shared::Triangulation<dim>::partition_zorder);

    test<dim>(triangulation, triangulation_other);
  }
  deallog.pop();
}
