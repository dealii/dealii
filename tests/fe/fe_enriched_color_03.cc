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
 * Test function - ColorEnriched::internal::set_cellwise_color_set_and_fe_index
 * for a set of predicates.
 * Check for each cell, if appropriate FE index and color-index map is set.
 * Color-index map associates different colors of different enrichment
 * functions with corresponding enrichment function index.
 */

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

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

  /*
   * Run through active cells to check FE index, colors of
   * enrichment functions associated with it.
   *
   * A unique color set corresponds to an FE index.
   *
   * Eg: If an FE index 1 corresponds to color set {2,3},
   * means that a cell with FE index 1 has enrichment functions
   * which are colored 2 and 3. Here different enrichment function
   * have same color 2 but for a given cell only one of them would
   * be relevant. So all additional information we need is which
   * enrichment function is relevant that has color 2. We need to
   * do the same thing with color 2 as well.
   *
   * Each cell is assigned unique material id by the function
   * set_cellwise_color_set_and_fe_index. Now using material id,
   * each cell is associated with a map which assigns a color to a
   * particular enrichment function id.
   */
  auto cell = dof_handler.begin_active();
  auto endc = dof_handler.end();
  for (unsigned int cell_index = 0; cell != endc; ++cell, ++cell_index)
    {
      // print true predicates for a cell as a binary code
      // 0 1 0 indicates predicate with index 2 is true in the cell.
      deallog << cell->index() << ":predicates=";
      for (auto predicate : vec_predicates)
        deallog << predicate(cell) << ':';

      // print color predicate pairs for a cell
      // Note that here material id is used to identify cells
      // Here (1,2) indicates predicate 2 of color 1 is relevant for cell.
      deallog << "(color, enrichment_func_id):";
      for (auto color_predicate_pair :
           cellwise_color_predicate_map[cell->material_id()])
        {
          deallog << '(' << color_predicate_pair.first << ','
                  << color_predicate_pair.second << "):";
        }

      // For a cell, print FE active index and corresponding FE set.
      //{1,2} indicates 2 enrichment functions of color 1 and 2 are relevant.
      deallog << ":fe_active_index:" << cell->active_fe_index() << ":fe_set:";
      for (auto fe_set_element : fe_sets[cell->active_fe_index()])
        deallog << fe_set_element << ':';
      deallog << std::endl;
    }

  return 0;
}
