// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a serial triangulation with periodic face in x-direction and copy it.

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

#include "../grid/tests.h"


template <int dim>
void
test(const int n_refinements, const int n_subdivisions, MPI_Comm comm)
{
  const double left  = 0;
  const double right = 1;

  auto add_periodicity = [&](dealii::Triangulation<dim> &tria) {
    std::vector<
      GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
         periodic_faces;
    auto cell = tria.begin();
    auto endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int face_number : GeometryInfo<dim>::face_indices())
        if (std::fabs(cell->face(face_number)->center()[0] - left) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(1);
        else if (std::fabs(cell->face(face_number)->center()[0] - right) <
                 1e-12)
          cell->face(face_number)->set_all_boundary_ids(2);
        else if (dim >= 2 &&
                 std::fabs(cell->face(face_number)->center()[1] - left) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(3);
        else if (dim >= 2 && std::fabs(cell->face(face_number)->center()[1] -
                                       right) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(4);
        else if (dim >= 3 &&
                 std::fabs(cell->face(face_number)->center()[2] - left) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(5);
        else if (dim >= 3 && std::fabs(cell->face(face_number)->center()[2] -
                                       right) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(6);

    if (dim >= 1)
      GridTools::collect_periodic_faces(tria, 1, 2, 0, periodic_faces);
    if (dim >= 2)
      GridTools::collect_periodic_faces(tria, 3, 4, 1, periodic_faces);
    if (dim >= 3)
      GridTools::collect_periodic_faces(tria, 5, 6, 2, periodic_faces);

    tria.add_periodicity(periodic_faces);
  };

  // create serial triangulation
  Triangulation<dim> basetria;
  GridGenerator::subdivided_hyper_cube(basetria, n_subdivisions);
  // new: add periodicity on serial mesh
  add_periodicity(basetria);
  basetria.refine_global(n_refinements);

  GridTools::partition_triangulation_zorder(
    Utilities::MPI::n_mpi_processes(comm), basetria);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  // new: add periodicity on fullydistributed mesh (!!!)
  add_periodicity(tria_pft);

  // test triangulation
  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(tria_pft);
  dof_handler.distribute_dofs(fe);

  // print statistics
  print_statistics(tria_pft);
  print_statistics(dof_handler);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  {
    deallog.push("1d");
    const int n_refinements  = 4;
    const int n_subdivisions = 8;
    test<1>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }

  {
    deallog.push("2d");
    const int n_refinements  = 3;
    const int n_subdivisions = 8;
    test<2>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }

  {
    deallog.push("3d");
    const int n_refinements  = 3;
    const int n_subdivisions = 4;
    test<3>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }
}
