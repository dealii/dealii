// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
//

// test DoFTools::extract_dofs_with_support_contained_within() for the
// p::d::Tria. As an MPI test, make sure that the total number of shape
// functions that are non-zero within the domain is the same.

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

#include <list>
#include <set>
#include <sstream>

#include "../tests.h"


template <int dim>
bool
pred_d(const typename DoFHandler<dim>::active_cell_iterator &cell)
{
  return (cell->center()[0] < 0.5 && cell->center()[1] < 0.5);
}


std::string
output_name(const unsigned int flag, const unsigned int subdomain)
{
  return "output" + Utilities::int_to_string(flag) + "_" +
         Utilities::int_to_string(subdomain) + ".vtu";
}


template <int dim>
void
test(const unsigned int flag)
{
  // Setup system
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);

  GridGenerator::hyper_rectangle(triangulation,
                                 Point<dim>(0, 0),
                                 Point<dim>(1, 1));

  if (flag == 0)
    triangulation.refine_global(2);
  else
    triangulation.refine_global(1);

  DoFHandler<dim> dh(triangulation);

  // Extra refinement to generate hanging nodes
  for (typename DoFHandler<dim>::active_cell_iterator cell = dh.begin_active();
       cell != dh.end();
       ++cell)
    if (cell->is_locally_owned() &&
        ((flag == 1 && pred_d<dim>(cell)) || (flag == 2 && !pred_d<dim>(cell))))
      cell->set_refine_flag();

  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement();

  if (flag > 0)
    triangulation.refine_global(1);

  FE_Q<dim> fe(2);
  dh.distribute_dofs(fe);

  const IndexSet locally_relevant_set =
    DoFTools::extract_locally_relevant_dofs(dh);

  AffineConstraints<double> cm;
  cm.reinit(dh.locally_owned_dofs(), locally_relevant_set);
  DoFTools::make_hanging_node_constraints(dh, cm);
  cm.close();

  const IndexSet support = DoFTools::extract_dofs_with_support_contained_within(
    dh,
    std::function<bool(const typename DoFHandler<dim>::active_cell_iterator &)>(
      &pred_d<dim>),
    cm);
  const IndexSet support_local = support & dh.locally_owned_dofs();

  deallog << support.n_elements() << std::endl;

  // now accumulate the number of indices and make sure it's the same for
  // various runs with different number of MPI cores
  const unsigned int dofs_support =
    Utilities::MPI::sum(support_local.n_elements(), MPI_COMM_WORLD);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      deallog << "accumulated: " << std::endl;
      deallog << dofs_support << std::endl;
    }

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
        output_name(flag, triangulation.locally_owned_subdomain());

      std::ofstream output(filename);
      data_out.write_vtu(output);

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          std::vector<std::string> filenames;
          for (unsigned int i = 0;
               i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
               ++i)
            filenames.push_back(output_name(flag, i));

          const std::string pvtu_filename =
            "output" + Utilities::int_to_string(flag) + ".pvtu";
          std::ofstream pvtu_output(pvtu_filename);
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

  test<2>(0);
  test<2>(1);
  test<2>(2);

  deallog.pop();
}
