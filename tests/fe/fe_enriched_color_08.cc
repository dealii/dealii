// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
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
 * Test if correct constraints are assigned for hp::DoFHandler using
 * hp::FECollection constructed by ColorEnriched::Helper.
 *
 * Due to the bug #1496 (https://github.com/dealii/dealii/issues/1496),
 * the constraints at certain interfaces may not be correct.
 * This tests if special treatment of the interface between enriched cells
 * results in the correct constraints despite the bug.
 */

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>

#include <map>

#include "../tests.h"

// used only for debugging
template <int dim>
void
plot_shape_function(DoFHandler<dim> &dof_handler, unsigned int patches = 5)
{
  deallog << "...start plotting shape function" << std::endl;
  deallog << "Patches for output: " << patches << std::endl;

  AffineConstraints<double> constraints;
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  // find set of dofs which belong to enriched cells
  std::set<unsigned int> enriched_cell_dofs;
  for (auto &cell : dof_handler.active_cell_iterators())
    if (cell->active_fe_index() != 0)
      {
        unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        enriched_cell_dofs.insert(local_dof_indices.begin(),
                                  local_dof_indices.end());
      }

  // output to check if all is good:
  std::vector<Vector<double>> shape_functions;
  std::vector<std::string>    names;
  for (auto dof : enriched_cell_dofs)
    {
      Vector<double> shape_function;
      shape_function.reinit(dof_handler.n_dofs());
      shape_function[dof] = 1.0;

      // if the dof is constrained, first output unconstrained vector
      names.push_back(std::string("C_") +
                      dealii::Utilities::int_to_string(dof, 2));
      shape_functions.push_back(shape_function);

      //      names.push_back(std::string("UC_") +
      //                      dealii::Utilities::int_to_string(s,2));

      //      // make continuous/constraint:
      //      constraints.distribute(shape_function);
      //      shape_functions.push_back(shape_function);
    }

  if (dof_handler.n_dofs() < 100)
    {
      deallog << "...start printing support points" << std::endl;

      std::map<types::global_dof_index, Point<dim>> support_points;
      MappingQ1<dim>                                mapping;
      hp::MappingCollection<dim>                    hp_mapping;
      for (unsigned int i = 0; i < dof_handler.get_fe_collection().size(); ++i)
        hp_mapping.push_back(mapping);
      DoFTools::map_dofs_to_support_points(hp_mapping,
                                           dof_handler,
                                           support_points);

      const std::string base_filename =
        "DOFs" + dealii::Utilities::int_to_string(dim) + "_p" +
        dealii::Utilities::int_to_string(0);

      const std::string filename = base_filename + ".gp";
      std::ofstream     f(filename);

      f << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
        << std::endl
        << "set output \"" << base_filename << ".png\"" << std::endl
        << "set size square" << std::endl
        << "set view equal xy" << std::endl
        << "unset xtics                                                          "
           "                         "
        << std::endl
        << "unset ytics" << std::endl
        << "unset grid" << std::endl
        << "unset border" << std::endl
        << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 "
           "offset 1,1 notitle"
        << std::endl;
      GridOut grid_out;
      grid_out.write_gnuplot(dof_handler.get_triangulation(), f);
      f << 'e' << std::endl;

      DoFTools::write_gnuplot_dof_support_point_info(f, support_points);

      f << 'e' << std::endl;

      deallog << "...finished printing support points" << std::endl;
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  // get material ids:
  Vector<float> fe_index(dof_handler.get_triangulation().n_active_cells());
  for (auto &cell : dof_handler.active_cell_iterators())
    {
      fe_index[cell->active_cell_index()] = cell->active_fe_index();
    }
  data_out.add_data_vector(fe_index, "fe_index");

  for (unsigned int i = 0; i < shape_functions.size(); ++i)
    data_out.add_data_vector(shape_functions[i], names[i]);

  data_out.build_patches(patches);

  std::string   filename = "shape_functions.vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);

  deallog << "...finished plotting shape functions" << std::endl;
}



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
    return ((i->center() - origin).norm_square() <= radius * radius);
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
  Point<dim>         p1(0, 0), p2(3, 1);
  GridGenerator::subdivided_hyper_rectangle(triangulation, {3, 1}, p1, p2);

  // Make predicates resulting in three adjacent domains
  // The relevant enrichment functions for the domains are [][1,3][2,3]
  // [] means no enrichment function.
  std::vector<predicate_function<dim>> vec_predicates;
  vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(1.5, 0.5), 0.5));
  vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(2.5, 0.5), 0.5));
  vec_predicates.push_back(EnrichmentPredicate<dim>(Point<dim>(2, 0.5), 1));

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
  FE_Q<dim> fe_base(2);
  FE_Q<dim> fe_enriched(1);

  static ColorEnriched::Helper<dim> fe_space(fe_base,
                                             fe_enriched,
                                             vec_predicates,
                                             vec_enrichments);
  hp::FECollection<dim>             fe_collection(
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

  deallog << "Face dominating set for 1 and 2: "
          << fe_collection.find_dominating_fe_extended({1, 2}, /*codim=*/1)
          << std::endl;

  dof_handler.distribute_dofs(fe_collection);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  constraints.print(deallog.get_file_stream());

  //  plot_shape_function(dof_handler);

  dof_handler.clear();
  return 0;
}
