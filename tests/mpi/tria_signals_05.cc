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

// This test is based on tria_signals_02. The difference is in this case we know
// that p4est is doing mesh smoothing beyond class Triangulation. The case setup
// is borrowed from tests/distributed_grids/2d_refinement_10.

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

  TriaType tria(MPI_COMM_WORLD);

  {
    std::vector<unsigned int> repetitions;
    Point<dim>                p1;
    Point<dim>                p2;

    for (unsigned int d = 0; d < dim; ++d)
      {
        repetitions.push_back(2);
        p1[d] = 0.0;
        p2[d] = 1.0;
      }
    GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);
  }

  SignalListener<dim, spacedim> count_cell_via_signal(tria);

  for (unsigned int n_loop = 1; n_loop < 5; ++n_loop)
    {
      {
        Point<dim> p;
        for (unsigned int d = 0; d < dim; ++d)
          {
            p[d] = 0.5 - std::pow(0.5, 1.0 + n_loop);
          }
        typename TriaType::active_cell_iterator cell = tria.begin_active();
        for (; cell != tria.end(); ++cell)
          if (cell->is_locally_owned() && ((cell->center()).distance(p) < 1e-4))
            {
              cell->set_refine_flag();
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

  // GridGenerator::subdivided_hyper_rectangle do not accept
  // parallel::distributed::Triangulation<2, 3>.

  {
    const int dim      = 3;
    const int spacedim = 3;
    test<dim, spacedim>();
  }

  return (0);
}
