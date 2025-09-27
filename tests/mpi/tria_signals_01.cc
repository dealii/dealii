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


// test signals of Triangulation class in parallel

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

int signal_counter_create                  = 0;
int signal_counter_pre_refinement          = 0;
int signal_counter_post_refinement         = 0;
int signal_counter_pre_coarsening_on_cell  = 0;
int signal_counter_post_refinement_on_cell = 0;
int signal_counter_copy                    = 0;
int signal_counter_clear                   = 0;
int signal_counter_any_change              = 0;

std::ofstream logfile;

void
f_create()
{
  ++signal_counter_create;
  return;
}

void
f_pre_refinement()
{
  ++signal_counter_pre_refinement;
  return;
}

void
f_post_refinement()
{
  ++signal_counter_post_refinement;
  return;
}

template <int dim, int spacedim>
void
f_pre_coarsening_on_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &)
{
  ++signal_counter_pre_coarsening_on_cell;
  return;
}

template <int dim, int spacedim>
void
f_post_refinement_on_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &)
{
  ++signal_counter_post_refinement_on_cell;
  return;
}

template <int dim, int spacedim>
void
f_copy(const Triangulation<dim, spacedim> &)
{
  ++signal_counter_copy;
  return;
}

void
f_clear()
{
  ++signal_counter_clear;
  return;
}

void
f_any_change()
{
  ++signal_counter_any_change;
  return;
}


template <int dim, int spacedim>
void
test()
{
  signal_counter_create                  = 0;
  signal_counter_pre_refinement          = 0;
  signal_counter_post_refinement         = 0;
  signal_counter_pre_coarsening_on_cell  = 0;
  signal_counter_post_refinement_on_cell = 0;
  signal_counter_copy                    = 0;
  signal_counter_clear                   = 0;
  signal_counter_any_change              = 0;

  using TriaType = parallel::distributed::Triangulation<dim, spacedim>;

  {
    const std::string prefix = Utilities::int_to_string(dim, 1) + "d-" +
                               Utilities::int_to_string(spacedim, 1) + "d";
    deallog.push(prefix.c_str());
  }

  TriaType tria(MPI_COMM_WORLD);

  tria.signals.create.connect(&f_create);
  tria.signals.pre_refinement.connect(&f_pre_refinement);
  tria.signals.post_refinement.connect(&f_post_refinement);
  tria.signals.pre_coarsening_on_cell.connect(
    &f_pre_coarsening_on_cell<dim, spacedim>);
  tria.signals.post_refinement_on_cell.connect(
    &f_post_refinement_on_cell<dim, spacedim>);
  tria.signals.copy.connect(&f_copy<dim, spacedim>);
  tria.signals.clear.connect(&f_clear);
  tria.signals.any_change.connect(&f_any_change);

  GridGenerator::hyper_cube(tria);

  // Test signal on copying
  {
    TriaType tria_cpoier(MPI_COMM_WORLD);
    tria_cpoier.copy_triangulation(tria);
  }

  tria.refine_global(2);

  // Test signal on coarsening
  {
    typename TriaType::active_cell_iterator       cell = tria.begin_active();
    const typename TriaType::active_cell_iterator endc = tria.end();

    for (; cell != endc; ++cell)
      {
        cell->set_coarsen_flag();
      }
    tria.execute_coarsening_and_refinement();
  }

  tria.clear();

  deallog << "n_signal_create : " << signal_counter_create << std::endl;
  deallog << "n_signal_pre_refinement : " << signal_counter_pre_refinement
          << std::endl;
  deallog << "n_signal_post_refinement : " << signal_counter_post_refinement
          << std::endl;
  deallog << "n_signal_pre_coarsening_on_cell : "
          << signal_counter_pre_coarsening_on_cell << std::endl;
  deallog << "n_signal_post_refinement_on_cell : "
          << signal_counter_post_refinement_on_cell << std::endl;
  deallog << "n_signal_copy : " << signal_counter_copy << std::endl;
  deallog << "n_signal_clear : " << signal_counter_clear << std::endl;
  deallog << "n_signal_any_change : " << signal_counter_any_change << std::endl;

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
