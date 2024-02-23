// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Test the function ColorEnriched::internal::color_predicates.
 * Check if the function correctly colors vector of predicates
 * using graph coloring.
 *
 * Two predicates are said to be connected if cells belonging to
 * different predicates touch each other.
 */

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_enriched.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

/*
 * Predicate function needed by ColorEnriched::internal::color_predicates
 * implemented using a struct.
 */
template <int dim>
struct EnrichmentPredicate
{
  EnrichmentPredicate(const Point<dim> origin, const double radius)
    : origin(origin)
    , radius(radius)
  {}

  template <class Iterator>
  bool
  operator()(const Iterator &i) const
  {
    return ((i->center() - origin).norm_square() < radius * radius);
  }

  const Point<dim> &
  get_origin()
  {
    return origin;
  }
  const double &
  get_radius()
  {
    return radius;
  }

private:
  const Point<dim> origin;
  const double     radius;
};



/*
 * Type used to defined vector of predicates needed by the function
 * ColorEnriched::internal::color_predicates.
 */
template <int dim>
using predicate_function =
  std::function<bool(const typename Triangulation<dim>::cell_iterator &)>;



int
main(int argc, char **argv)
{
  // Initialize MPI as required by Zoltan library used for graph coloring by
  // this test.
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  MPILogInitAll                    all;

  // Make basic grid
  const unsigned int dim = 2;
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -20, 20);
  triangulation.refine_global(4);
  DoFHandler<dim> dof_handler(triangulation);

  // check the coloring function on different set of predicates.
  std::vector<predicate_function<dim>> vec_predicates;
  std::vector<unsigned int>            predicate_colors;
  {
    // case 1: predicates are not connected
    vec_predicates.clear();
    vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(-10, 10), 2));
    vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(0, 0), 2));

    predicate_colors.resize(vec_predicates.size());

    ColorEnriched::internal::color_predicates<dim>(dof_handler,
                                                   vec_predicates,
                                                   predicate_colors);

    deallog << "Case 1" << std::endl;
    for (auto i : predicate_colors)
      {
        deallog << i << std::endl;
      }
  }

  {
    // case 2: Two predicates that are connected.
    vec_predicates.clear();
    vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(-10, 10), 2));
    vec_predicates.push_back(
      EnrichmentPredicate<dim>(Point<dim>(-7.5, 7.5), 2));

    predicate_colors.resize(vec_predicates.size());

    ColorEnriched::internal::color_predicates(dof_handler,
                                              vec_predicates,
                                              predicate_colors);

    deallog << "Case 2" << std::endl;
    for (auto i : predicate_colors)
      {
        deallog << i << std::endl;
      }
  }

  {
    // case 3: connections between predicates is as follows -
    // 0-1 (overlap connection),
    // 3-4 (edge connection)
    vec_predicates.clear();
    vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(-10, 10), 2));
    vec_predicates.push_back(
      EnrichmentPredicate<dim>(Point<dim>(-7.5, 7.5), 2));
    vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(0, 0), 2));
    vec_predicates.push_back(
      EnrichmentPredicate<dim>(Point<dim>(7.5, -7.5), 2));
    vec_predicates.push_back(
      EnrichmentPredicate<dim>(Point<dim>(12.5, -12.5), 2));

    predicate_colors.resize(vec_predicates.size());

    ColorEnriched::internal::color_predicates(dof_handler,
                                              vec_predicates,
                                              predicate_colors);

    deallog << "Case 3" << std::endl;
    for (auto i : predicate_colors)
      {
        deallog << i << std::endl;
      }
  }
  return 0;
}
