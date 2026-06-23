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
//

// test DoFTools::extract_dofs_with_support_contained_within() for the
// p::d::Tria. As an MPI test, make sure that the total number of shape
// functions that are non-zero within the predicate domain is the same
// regardless of how many MPI processes are used.
// Note: the function DoFTools::extract_dofs_with_support_contained_within()
// does not return the total number of shape functions contained within the
// predicate domain, but only a local subset that contains locally active DoFs,
// as well as DoFs by which a locally active DoF is constrained.
/*
    This test case considers the special case of a locally active Dof (DoF
  marked with an x) that constrains a DoF (DoF marked with an o), which is
  supported on predicate ghost cells (marked by #g#) of the current process and
  on a artificial cell the current process has no information about.

    +---+---+---+---+
    |###|###|#g#| a1|
    +---+---+---+---+
    |###|###|#g#| a2|
    +---+---x---o---+
    |#######|###g###|
    |#######|###g###|
    +-------+-------+
    current |neighboring
    process |process

            ##### = predicate region
            #g#   = predicate ghost cell
          a1,a2  = artificial cell
*/

#include "deal.II/base/index_set.h"
#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

void
test()
{
  constexpr unsigned                        dim = 2;
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            {2, 1},
                                            Point<dim>(0, 0),
                                            Point<dim>(2.0, 1.0));

  triangulation.refine_global(1);

  // Create a triangulation that partitions right down the middle, where the
  // upper half is refined twice and the upper quarter is refined three times
  // resulting in hanging node constraints within the ghost layer.

  const auto predicate_upper_half =
    [](const typename Triangulation<dim>::active_cell_iterator &cell) {
      return cell->center()[1] > 0.5;
    };

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned() && predicate_upper_half(cell))
        {
          cell->set_refine_flag();
        }
    }
  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement();

  const auto predicate_upper_quarter =
    [](const typename Triangulation<dim>::active_cell_iterator &cell) {
      return cell->center()[1] > 0.75;
    };
  for (const auto &cell : triangulation.active_cell_iterators())

    {
      if (cell->is_locally_owned() && predicate_upper_quarter(cell))
        {
          cell->set_refine_flag();
        }
    }

  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement();



  DoFHandler<dim> dh(triangulation);

  FE_Q<dim> fe(1);
  dh.distribute_dofs(fe);

  const IndexSet &locally_owned_set = dh.locally_owned_dofs();
  const IndexSet  locally_relevant_set =
    DoFTools::extract_locally_relevant_dofs(dh);

  AffineConstraints<double> cm;
  cm.reinit(locally_owned_set, locally_relevant_set);
  DoFTools::make_hanging_node_constraints(dh, cm);
  cm.close();

  const IndexSet locally_active = DoFTools::extract_locally_active_dofs(dh);

  deallog << "Number of locally active dofs: " << locally_active.n_elements()
          << std::endl;

  // Choose a predicate such that we end up with a locally active DoF
  // constraining a ghost DoF, where the ghost DoF is supported outside the
  // known domain (locally owned cells and ghost cells) of the MPI process.
  const auto support_predicate =
    [&predicate_upper_half, &predicate_upper_quarter](
      const typename Triangulation<dim>::active_cell_iterator &cell) {
      return (cell->center()[0] > 0.5 && cell->center()[0] < 1.25 &&
              predicate_upper_half(cell) && cell->center()[1] < 0.75) ||
             (cell->center()[0] > 0.675 && cell->center()[0] < 1.125 &&
              predicate_upper_quarter(cell) && cell->center()[1] < 1.0);
    };

  const IndexSet support =
    DoFTools::extract_dofs_with_support_contained_within(dh,
                                                         support_predicate,
                                                         cm);

  const IndexSet support_locally_owned = support & dh.locally_owned_dofs();

  deallog << "Number of dofs with contained support: " << support.n_elements()
          << std::endl;
  deallog << "Number of owned dofs with contained support: "
          << support_locally_owned.n_elements() << std::endl;

  const unsigned int dofs_support =
    Utilities::MPI::sum(support_locally_owned.n_elements(), MPI_COMM_WORLD);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      deallog << "Total number of dofs with contained support: " << dofs_support
              << std::endl;
    }
  deallog << std::endl;

  // print grid and DoFs for visual inspection
  if (false)
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler(dh);

      std::vector<LinearAlgebra::distributed::Vector<double>> shape_functions(
        dh.n_dofs());
      for (unsigned int i = 0; i < dh.n_dofs(); ++i)
        {
          LinearAlgebra::distributed::Vector<double> &s = shape_functions[i];
          s.reinit(dh.locally_owned_dofs(),
                   locally_relevant_set,
                   MPI_COMM_WORLD);
          s = 0.;
          if (dh.locally_owned_dofs().is_element(i))
            s[i] = 1.0;
          s.compress(VectorOperation::insert);
          cm.distribute(s);
          s.update_ghost_values();

          data_out.add_data_vector(s,
                                   std::string("N_") +
                                     Utilities::int_to_string(i));
        }

      Vector<float> subdomain(triangulation.n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = triangulation.locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");
      data_out.build_patches();

      const std::string filename =
        "output_subdomain_" +
        std::to_string(triangulation.locally_owned_subdomain()) + ".vtu";

      std::ofstream output(filename);
      data_out.write_vtu(output);

      const unsigned rank    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      const unsigned n_ranks = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      if (rank == 0)
        {
          std::vector<std::string> filenames;
          filenames.reserve(n_ranks);
          for (unsigned int rank = 0; rank < n_ranks; ++rank)
            filenames.push_back("output_subdomain_" + std::to_string(rank) +
                                ".vtu");

          const std::string pvtu_filename = "output.pvtu";
          std::ofstream     pvtu_output(pvtu_filename);
          data_out.write_pvtu_record(pvtu_output, filenames);
        }

      for (unsigned int i = 0; i < n_ranks; ++i)
        {
          MPI_Barrier(MPI_COMM_WORLD);
          if (i == rank)
            {
              std::cout << "-------------------- " << rank << std::endl;
              cm.print(std::cout);

              std::cout << "local support:" << std::endl;
              (support & locally_owned_set).print(std::cout);
              std::cout << "support:" << std::endl;
              support.print(std::cout);
              std::cout << "locally owned:" << std::endl;
              locally_owned_set.print(std::cout);
              std::cout << "locally relevant set:" << std::endl;
              locally_relevant_set.print(std::cout);
            }
        }

      std::map<types::global_dof_index, Point<dim>> support_points;
      MappingQ1<dim>                                mapping;
      DoFTools::map_dofs_to_support_points(mapping, dh, support_points);

      const std::string filename_gnuplot = "grid_n_ranks_" +
                                           std::to_string(n_ranks) + "_rank_" +
                                           std::to_string(rank);
      std::ofstream f(filename_gnuplot + ".gp");

      f << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
        << std::endl
        << "set output \"" << filename_gnuplot << ".png\"" << std::endl
        << "set size square" << std::endl
        << "set view equal xy" << std::endl
        << "unset xtics" << std::endl
        << "unset ytics" << std::endl
        << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 1,1 notitle"
        << std::endl;
      GridOut().write_gnuplot(triangulation, f);
      f << 'e' << std::endl;

      DoFTools::write_gnuplot_dof_support_point_info(f, support_points);
      f << 'e' << std::endl;
    }

  dh.clear();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  MPILogInitAll log;

  deallog.push(Utilities::int_to_string(myid));

  test();

  deallog.pop();
}
