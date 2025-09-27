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

// uncomment when debugging
// #define DATA_OUT_FE_ENRICHED

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



template <int dim>
void
plot_shape_function(DoFHandler<dim> &dof_handler, unsigned int patches = 5)
{
  std::cout << "n_cells: " << dof_handler.get_triangulation().n_active_cells()
            << std::endl;

  AffineConstraints<double> constraints;
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  // output to check if all is good:
  std::vector<Vector<double>> shape_functions;
  std::vector<std::string>    names;
  for (unsigned int s = 0; s < dof_handler.n_dofs(); ++s)
    {
      Vector<double> shape_function;
      shape_function.reinit(dof_handler.n_dofs());
      shape_function[s] = 1.0;

      // if the dof is constrained, first output unconstrained vector
      if (constraints.is_constrained(s))
        {
          names.push_back(std::string("UN_") +
                          dealii::Utilities::int_to_string(s, 2));
          shape_functions.push_back(shape_function);
        }

      names.push_back(std::string("N_") +
                      dealii::Utilities::int_to_string(s, 2));

      // make continuous/constrain:
      constraints.distribute(shape_function);
      shape_functions.push_back(shape_function);
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  // get material ids:
  Vector<float> fe_index(dof_handler.get_triangulation().n_active_cells());
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (unsigned int index = 0; cell != endc; ++cell, ++index)
    {
      fe_index[index] = cell->active_fe_index();
    }
  data_out.add_data_vector(fe_index, "fe_index");

  for (unsigned int i = 0; i < shape_functions.size(); ++i)
    data_out.add_data_vector(shape_functions[i], names[i]);

  data_out.build_patches(patches);

  std::string filename =
    "hp-shape_functions_" + dealii::Utilities::int_to_string(dim) + "D.vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);
}



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

  // find colors for predicates
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
      // constant function.
      Functions::ConstantFunction<dim> func(10 + i); // constant function
      vec_enrichments.push_back(
        std::make_shared<Functions::ConstantFunction<dim>>(func));
    }

  // Construct container for color enrichment functions needed
  // by function make_colorwise_enrichment_functions
  std::vector<std::function<const Function<dim> *(
    const typename Triangulation<dim>::cell_iterator &)>>
    color_enrichments;

  ColorEnriched::internal::make_colorwise_enrichment_functions<dim, dim>(
    num_colors,      // needs number of colors
    vec_enrichments, // enrichment functions based on predicate id
    cellwise_color_predicate_map,
    color_enrichments);

  // Construct object needed to call make_fe_collection_from_colored_enrichments
  FE_Q<dim>             fe_base(2);
  FE_Q<dim>             fe_enriched(1);
  FE_Nothing<dim>       fe_nothing(1, true);
  hp::FECollection<dim> fe_collection;
  ColorEnriched::internal::make_fe_collection_from_colored_enrichments(
    num_colors,
    fe_sets,           // total list of color sets possible
    color_enrichments, // color wise enrichment functions
    fe_base,           // basic fe element
    fe_enriched,       // fe element multiplied by enrichment function
    fe_nothing,
    fe_collection);

  // print all the different FE sets needed by different cells
  deallog << "fe sets:" << std::endl;
  for (auto fe_set : fe_sets)
    {
      deallog << "color:";
      for (auto color : fe_set)
        deallog << ':' << color;
      deallog << std::endl;
    }

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

#ifdef DATA_OUT_FE_ENRICHED
  GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(
                                       MPI_COMM_WORLD),
                                     triangulation,
                                     SparsityTools::Partitioner::zoltan);
  dof_handler.distribute_dofs(*fe_collection);

  plot_shape_function<dim>(dof_handler, 5);
#endif

  dof_handler.clear();
  return 0;
}
