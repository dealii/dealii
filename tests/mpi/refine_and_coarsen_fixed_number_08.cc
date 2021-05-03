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


// verify the binary search algorithm in different compositions of
// criteria and refinement and coarsening fractions for
// parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number().

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
verify(parallel::distributed::Triangulation<dim> &tr,
       const Vector<float> &                      criteria,
       const float                                refinement_fraction,
       const float                                coarsening_fraction)
{
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
    tr, criteria, refinement_fraction, coarsening_fraction);

  unsigned int n_refine_flags = 0, n_coarsen_flags = 0;
  for (const auto &cell : tr.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        if (cell->refine_flag_set())
          {
            ++n_refine_flags;
            cell->clear_refine_flag();
          }
        if (cell->coarsen_flag_set())
          {
            ++n_coarsen_flags;
            cell->clear_coarsen_flag();
          }
      }

  const unsigned int n_global_refine_flags =
                       Utilities::MPI::sum(n_refine_flags, MPI_COMM_WORLD),
                     n_global_coarsen_flags =
                       Utilities::MPI::sum(n_coarsen_flags, MPI_COMM_WORLD);

  deallog << "  refinement_fraction:" << refinement_fraction
          << " coarsening_fraction:" << coarsening_fraction << std::endl
          << "    n_refine_flags:" << n_global_refine_flags
          << " n_coarsen_flags:" << n_global_coarsen_flags << std::endl;
}


template <int dim>
void
test()
{
  const unsigned int n_cells             = 100;
  const float        refinement_fraction = 0.1, coarsening_fraction = 0.1;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  std::vector<unsigned int>                 rep(dim, 1);
  rep[0] = n_cells;
  Point<dim> p1, p2;
  for (unsigned int d = 0; d < dim; ++d)
    {
      p1[d] = 0;
      p2[d] = (d == 0) ? n_cells : 1;
    }
  GridGenerator::subdivided_hyper_rectangle(tr, rep, p1, p2);

  deallog << "n_cells:" << n_cells << std::endl;
  Vector<float> criteria(n_cells);
  {
    deallog << "criteria:[b>0,e>0]" << std::endl;

    std::iota(criteria.begin(), criteria.end(), 1);
    Assert(criteria(0) > 0. && criteria(n_cells - 1) > 0., ExcInternalError());

    verify<dim>(tr, criteria, refinement_fraction, coarsening_fraction);
    verify<dim>(tr, criteria, refinement_fraction, 0.);
    verify<dim>(tr, criteria, 0., coarsening_fraction);
  }

  {
    deallog << "criteria:[b=0,e>0]" << std::endl;

    std::iota(criteria.begin(), criteria.end(), 0);
    Assert(criteria(0) == 0. && criteria(n_cells - 1) > 0., ExcInternalError());

    verify<dim>(tr, criteria, refinement_fraction, coarsening_fraction);
    verify<dim>(tr, criteria, refinement_fraction, 0.);
    verify<dim>(tr, criteria, 0., coarsening_fraction);
  }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();
      deallog << std::setprecision(1);
    }

  test<2>();
}
