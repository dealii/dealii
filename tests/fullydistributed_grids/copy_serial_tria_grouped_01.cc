// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// Create the TriangulationDescription::Description with
// create_description_from_triangulation_in_groups, i.e. by a set of
// master processes.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <boost/archive/text_oarchive.hpp>

#include "./tests.h"

using namespace dealii;

template <int dim, int spacedim = dim>
void
test(int n_refinements, MPI_Comm comm)
{
  // 1) create TriangulationDescription::Description with
  // create_description_from_triangulation
  Triangulation<dim> basetria(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_L(basetria);
  basetria.refine_global(n_refinements);

  GridTools::partition_triangulation_zorder(
    Utilities::MPI::n_mpi_processes(comm), basetria);
  GridTools::partition_multigrid_levels(basetria);

  auto construction_data_1 =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm, true);

  // 2) create TriangulationDescription::Description with
  // create_description_from_triangulation_in_groups
  auto construction_data_2 = TriangulationDescription::Utilities::
    create_description_from_triangulation_in_groups<dim, spacedim>(
      [n_refinements](dealii::Triangulation<dim, spacedim> &tria) {
        GridGenerator::hyper_L(tria);
        tria.refine_global(n_refinements);
      },
      [](dealii::Triangulation<dim, spacedim> &tria,
         const MPI_Comm                        comm,
         const unsigned int /*group_size*/) {
        GridTools::partition_triangulation_zorder(
          Utilities::MPI::n_mpi_processes(comm), tria);
      },
      comm,
      3 /* group size */,
      true);

  // 3a) serialize first TriangulationDescription::Description and print
  {
    std::ostringstream            oss;
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << construction_data_1;
    deallog << oss.str() << std::endl;
  }

  // 3b) serialize second TriangulationDescription::Description and print
  {
    std::ostringstream            oss;
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << construction_data_2;
    deallog << oss.str() << std::endl;
  }

  // 4) the result has to be identical
  AssertThrow((construction_data_1 == construction_data_2 &&
               construction_data_1.comm == construction_data_2.comm),
              ExcMessage(
                "TriangulationDescription::Descriptions are not the same!"));
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const int      n_refinements = 3;
  const MPI_Comm comm          = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    test<2>(n_refinements, comm);
    deallog.pop();
  }
  {
    deallog.push("3d");
    test<3>(n_refinements, comm);
    deallog.pop();
  }
}
