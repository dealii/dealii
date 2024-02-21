// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <functional>
#include <ostream>

#include "../tests.h"

// Test on whether signals post_refinement_on_cell and pre_coarsening_on_cell
// could catch all cell changes.
// The test is designed to count cell number increase and decrease in signal
// calls and then compare the result against n_active_cells reported by Tria
// object. Absolute value change in n_active_cells is not concerned in this
// test.

// This test is based on tria_signals_03. The difference is here we have all
// smoothing flags enabled. Also increase refine fraction to one third to
// prevent refine flags getting smoothed out.

template <int dim, int spacedim>
class SignalListener
{
public:
  SignalListener(Triangulation<dim, spacedim> &tria_in)
    : n_active_cells(tria_in.n_active_cells())
    , tria(tria_in)
  {
    tria_in.signals.post_refinement_on_cell.connect(
      std::bind(&SignalListener<dim, spacedim>::count_on_refine,
                this,
                std::placeholders::_1));

    tria_in.signals.pre_coarsening_on_cell.connect(
      std::bind(&SignalListener<dim, spacedim>::count_on_coarsen,
                this,
                std::placeholders::_1));
  }

  int
  n_active_cell_gap()
  {
    return (n_active_cells - static_cast<int>(tria.n_active_cells()));
  }

private:
  void
  count_on_refine(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
  {
    n_active_cells += cell->n_children();
    --n_active_cells;

    return;
  }

  void
  count_on_coarsen(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
  {
    ++n_active_cells;
    n_active_cells -= cell->n_children();

    return;
  }

  int                                 n_active_cells;
  const Triangulation<dim, spacedim> &tria;
};


template <int dim, int spacedim>
void
test()
{
  using TriaType = parallel::distributed::Triangulation<dim, spacedim>;

  {
    const std::string prefix = Utilities::int_to_string(dim, 1) + "d-" +
                               Utilities::int_to_string(spacedim, 1) + "d";
    deallog.push(prefix.c_str());
  }

  // Option dealii::Triangulation<dim, spacedim>::maximum_smoothing can't
  // run in parallel at the time that this test is created.
  TriaType tria(
    MPI_COMM_WORLD,
    typename dealii::Triangulation<dim, spacedim>::MeshSmoothing(
      dealii::Triangulation<dim, spacedim>::smoothing_on_refinement |
      dealii::Triangulation<dim, spacedim>::smoothing_on_coarsening));


  GridGenerator::hyper_cube(tria);
  SignalListener<dim, spacedim> count_cell_via_signal(tria);

  tria.refine_global(1);


  // The following loop is borrowed from p4est_3d_refine_01 with some
  // modifications.
  for (int n_loop = 0;
       // Terminate loop on global information to prevent premature termination
       // on only part of processors. (n_loop < 20) is just a passive safety to
       // avoid infinite loop.
       (tria.n_global_active_cells() < 20000) && (n_loop < 20);
       ++n_loop)
    {
      std::vector<bool> flags(tria.n_active_cells(), false);

      // Refine one fifth of all cells each time (but at least one).
      // Note that only the own marked cells will be refined.
      // But refine flags on own cells could be effected by flags on ghost cells
      // through mesh smoothing.
      for (unsigned int i = 0; i < tria.n_active_cells() / 3 + 1; ++i)
        {
          const unsigned int x = Testing::rand() % flags.size();
          flags[x]             = true;
        }

      unsigned int index  = 0;
      unsigned int locals = 0;

      for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
             tria.begin_active();
           cell != tria.end();
           ++cell, ++index)
        if (flags[index])
          {
            if (cell->is_locally_owned())
              ++locals;
            cell->set_refine_flag();
          }

      if (locals > 5)
        {
          // Coarsen some cells randomly only if we have enough local cells
          // marked to be refined
          std::fill(flags.begin(), flags.end(), false);
          for (unsigned int i = 0; i < tria.n_active_cells() / 3; ++i)
            {
              const unsigned int x = Testing::rand() % flags.size();
              flags[x]             = true;
            }

          index = 0;
          for (typename Triangulation<dim, spacedim>::active_cell_iterator
                 cell = tria.begin_active();
               cell != tria.end();
               ++cell, ++index)
            if (flags[index] && !cell->refine_flag_set())
              {
                cell->set_coarsen_flag();
              }
        }

      tria.execute_coarsening_and_refinement();

      deallog << "n_loop: " << n_loop
              << ", n_cell_gap: " << count_cell_via_signal.n_active_cell_gap()
              << std::endl;
    }

  deallog.pop();
  return;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, /* int max_num_threads */ 1);
  MPILogInitAll log;

  // parallel::distributed::Triangulation<1, spacedim> is not valid.
  {
    const int dim      = 2;
    const int spacedim = 2;
    test<dim, spacedim>();
  }

  {
    const int dim      = 2;
    const int spacedim = 3;
    test<dim, spacedim>();
  }

  {
    const int dim      = 3;
    const int spacedim = 3;
    test<dim, spacedim>();
  }

  return (0);
}
