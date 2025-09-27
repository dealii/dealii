// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
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
 * Test function: ColorEnriched::internal
 * ::make_fe_collection_from_colored_enrichments for a set of predicates.
 *
 * The function return FE_Collection which is then printed to test.
 */

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>

#include <map>

#include "../tests.h"

/*
 * Testing helper class
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
  DoFHandler<dim>    dof_handler(triangulation);
  GridGenerator::hyper_cube(triangulation, -2, 2);
  triangulation.refine_global(2);

  // Make predicates. Predicate 0 and 1 overlap.
  std::vector<predicate_function<dim>> vec_predicates;
  vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(-1, 1), 1));
  vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(0, 1), 1));

  // Construct vector of enrichment functions
  std::vector<std::shared_ptr<Function<dim>>> vec_enrichments;
  vec_enrichments.reserve(vec_predicates.size());
  for (unsigned int i = 0; i < vec_predicates.size(); ++i)
    {
      // constant function.
      Functions::ConstantFunction<dim> func(10 + i); // constant function
      vec_enrichments.push_back(
        std::make_shared<Functions::ConstantFunction<dim>>(func));
    }

  // Construct helper class to construct FE collection
  FE_Q<dim>                         fe_base(2);
  FE_Q<dim>                         fe_enriched(1);
  static ColorEnriched::Helper<dim> fe_space(fe_base,
                                             fe_enriched,
                                             vec_predicates,
                                             vec_enrichments);
  const hp::FECollection<dim>      &fe_collection(
    fe_space.build_fe_collection(dof_handler));

  // check if fe_collection is correctly constructed by function
  deallog << "fe_collection[index] mapping:" << std::endl;
  for (unsigned int index = 0; index != fe_collection.size(); ++index)
    {
      deallog << "name:" << fe_collection[index].get_name() << std::endl;
      deallog << "n_blocks:" << fe_collection[index].n_blocks() << std::endl;
      deallog << "n_comp:" << fe_collection[index].n_components() << std::endl;
      deallog << "n_dofs:" << fe_collection[index].n_dofs_per_cell()
              << std::endl;
    }

  dof_handler.clear();
  return 0;
}
