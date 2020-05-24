// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// hp::Refinement::choose_p_over_h called future_fe_index_set() on cells that
// are not locally owned and triggered an assertion at some point.


#include <deal.II/distributed/shared_tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // setup
  const unsigned int n_cells = 2;

  parallel::shared::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_custom_signal);
  tr.signals.post_refinement.connect([&tr]() {
    // partition the triangulation by hand
    for (const auto &cell : tr.active_cell_iterators())
      cell->set_subdomain_id(cell->active_cell_index() %
                             Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD));
  });

  std::vector<unsigned int> rep(dim, 1);
  rep[0] = n_cells;
  Point<dim> p1, p2;
  for (unsigned int d = 0; d < dim; ++d)
    {
      p1[d] = 0;
      p2[d] = (d == 0) ? n_cells : 1;
    }
  GridGenerator::subdivided_hyper_rectangle(tr, rep, p1, p2);
  tr.refine_global(1);

  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 2; ++d)
    fes.push_back(FE_Q<dim>(d));

  hp::DoFHandler<dim> dh(tr);
  dh.set_fe(fes);

  // set flags
  for (auto cell = dh.begin(0); cell != dh.end(0); ++cell)
    {
      if (cell->id().to_string() == "0_0:")
        {
          // all children will be flagged for both h- and p-adaptation.
          // choose_p_over_h() will decide in favor of p-adaptation.
          for (unsigned int i = 0; i < cell->n_children(); ++i)
            {
              const auto &child = cell->child(i);
              if (child->is_locally_owned())
                {
                  child->set_future_fe_index(1);
                  child->set_coarsen_flag();
                }
            }
        }
      else if (cell->id().to_string() == "1_0:")
        {
          // all children will be flagged for both h-adaptation
          // and only one of them for p-adaptation.
          // choose_p_over_h() will decide in favor of h-adaptation.
          for (unsigned int i = 0; i < cell->n_children(); ++i)
            {
              const auto &child = cell->child(i);
              if (child->is_locally_owned())
                {
                  if (i == 0)
                    child->set_future_fe_index(1);
                  child->set_coarsen_flag();
                }
            }
        }
    }

  // decide between p and h flags
  hp::Refinement::choose_p_over_h(dh);

  // verify
  unsigned int h_flagged_cells = 0, p_flagged_cells = 0;
  for (const auto &cell : dh.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        if (cell->coarsen_flag_set())
          ++h_flagged_cells;
        if (cell->future_fe_index_set())
          ++p_flagged_cells;
      }
  const unsigned int global_h_flagged_cells =
                       Utilities::MPI::sum(h_flagged_cells, MPI_COMM_WORLD),
                     global_p_flagged_cells =
                       Utilities::MPI::sum(p_flagged_cells, MPI_COMM_WORLD);

  deallog << "h-flags:" << global_h_flagged_cells << std::endl
          << "p-flags:" << global_p_flagged_cells << std::endl
          << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }
}
