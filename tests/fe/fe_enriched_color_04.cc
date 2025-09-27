// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
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
 * Test function: ColorEnriched::internal::make_colorwise_enrichment_functions
 * for a set of predicates.
 *
 * The container of enrichment functions created by the function is
 * called color_enrichments because they return the correct enrichment
 * function if cell and color of the function are provided.
 */

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_enriched.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <map>

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
  DoFHandler<dim>    dof_handler(triangulation);
  GridGenerator::hyper_cube(triangulation, -2, 2);
  triangulation.refine_global(2);


  // Make predicates. Predicate 0 and 1 overlap.
  // Predicate 2 is connected to 0.
  std::vector<predicate_function<dim>> vec_predicates;
  vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(0, 1), 1));
  vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(-1, 1), 1));
  vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(1.5, -1.5), 1));

  // Do manual coloring since we are not testing coloring function here!
  std::vector<unsigned int> predicate_colors;
  predicate_colors.resize(vec_predicates.size());
  unsigned int num_colors =
    ColorEnriched::internal::color_predicates(dof_handler,
                                              vec_predicates,
                                              predicate_colors);


  // Make required objects to call function set_cellwise_color_set_and_fe_index
  std::map<unsigned int, std::map<unsigned int, unsigned int>>
                                      cellwise_color_predicate_map;
  std::vector<std::set<unsigned int>> fe_sets;
  ColorEnriched::internal::set_cellwise_color_set_and_fe_index(
    dof_handler,
    vec_predicates,
    predicate_colors,
    cellwise_color_predicate_map,
    fe_sets);

  // Construct vector of enrichment functions
  std::vector<std::shared_ptr<Function<dim>>> vec_enrichments;
  vec_enrichments.reserve(vec_predicates.size());
  for (unsigned int i = 0; i < vec_predicates.size(); ++i)
    {
      // constant function is chosen as enrichment function
      Functions::ConstantFunction<dim> func(i);
      vec_enrichments.push_back(
        std::make_shared<Functions::ConstantFunction<dim>>(func));
    }

  // Construct container for color enrichment functions needed
  // by function make_colorwise_enrichment_functions
  std::vector<std::function<const Function<dim> *(
    const typename Triangulation<dim>::cell_iterator &)>>
    color_enrichments;

  ColorEnriched::internal::make_colorwise_enrichment_functions<dim, dim>(
    num_colors,
    vec_enrichments,
    cellwise_color_predicate_map,
    color_enrichments);


  deallog << "color wise enrichment functions:" << std::endl;
  auto cell = dof_handler.begin_active();
  auto endc = dof_handler.end();
  for (unsigned int cell_index = 0; cell != endc; ++cell, ++cell_index)
    {
      // Print ids of predicates active in the cell
      unsigned int cell_id = cell->index();
      deallog << cell_id << ":predicates=";
      for (auto predicate : vec_predicates)
        deallog << predicate(cell) << ':';

      /*
       * Check if a color and enrichment index map exists for the cell.
       * If so print the color and corresponding enrichment functions
       * value at the cell center. Note that indices of color enrichment
       * functions starts with zero not one unlike color indices.
       */
      if (cellwise_color_predicate_map.count(cell->material_id()) == 1)
        for (unsigned int color = 1; color <= num_colors; ++color)
          if (cellwise_color_predicate_map.at(cell->material_id())
                .count(color) == 1)
            deallog << ":color:" << color << ":func_value:"
                    << color_enrichments[color - 1](cell)->value(
                         cell->center());

      deallog << std::endl;
    }
  return 0;
}
