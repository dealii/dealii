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

#include <deal.II/base/exception_macros.h>
#include <deal.II/base/index_set.h>
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

bool
predicate_left(const typename Triangulation<2, 2>::active_cell_iterator &cell)
{
  return (cell->center()[0] < 1.00);
}

bool
predicate_right(const typename Triangulation<2, 2>::active_cell_iterator &cell)
{
  return (cell->center()[0] > 1.00);
}

enum class Refine
{
  None,
  RightCell,
  LeftCell
};

std::string
to_string(Refine cell_to_refine)
{
  switch (cell_to_refine)
    {
      case Refine::None:
        return "none";
      case Refine::RightCell:
        return "right_cell";
      case Refine::LeftCell:
        return "left_cell";
    }
  DEAL_II_ASSERT_UNREACHABLE();
};


void
test(Refine cell_to_refine)
{
  constexpr unsigned                        dim = 2;
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            {2, 1},
                                            Point<dim>(0, 0),
                                            Point<dim>(2, 1));
  switch (cell_to_refine)
    {
      case Refine::None:
        break;
      case Refine::RightCell:
        for (const auto &cell : triangulation.active_cell_iterators())
          {
            if (cell->is_locally_owned() && predicate_right(cell))
              {
                cell->set_refine_flag();
              }
          }
        triangulation.signals.weight.connect(
          [](const typename Triangulation<dim>::cell_iterator &cell,
             const CellStatus                                  status) {
            // Give left unrefined cell a larger weight such that the mesh
            // partitions right down the middle. Consequently, process 0 owns
            // the unrefined cell on the left and process 1 the four refined
            // cells on the right.
            if (predicate_left(cell))
              return 4u;

            return 1u;
          });
        break;
      case Refine::LeftCell:
        for (const auto &cell : triangulation.active_cell_iterators())
          {
            if (cell->is_locally_owned() && predicate_left(cell))
              {
                cell->set_refine_flag();
              }
          }

        triangulation.signals.weight.connect(
          [](const typename Triangulation<dim>::cell_iterator &cell,
             const CellStatus                                  status) {
            // Give right unrefined cell a larger weight such that the mesh
            // partitions right down the middle. Consequently, process 0 owns
            // the four refined cells on the left and process 1 the one
            // unrefined cell on the right.
            if (predicate_right(cell))
              return 4u;

            return 1u;
          });
        break;
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

  deallog << "Refined cell: " << to_string(cell_to_refine) << std::endl;
  deallog << "Number of locally active dofs: " << locally_active.n_elements()
          << std::endl;

  const IndexSet support =
    DoFTools::extract_dofs_with_support_contained_within(dh,
                                                         &predicate_right,
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

      const std::string file_base_name =
        "output_refine_" + to_string(cell_to_refine);

      const std::string filename =
        file_base_name + "_subdomain_" +
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
            filenames.push_back(file_base_name + "_subdomain_" +
                                std::to_string(rank) + ".vtu");

          const std::string pvtu_filename = file_base_name + ".pvtu";
          std::ofstream     pvtu_output(pvtu_filename);
          data_out.write_pvtu_record(pvtu_output, filenames);
        }
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

  test(Refine::None);
  test(Refine::LeftCell);
  test(Refine::RightCell);

  deallog.pop();
}
