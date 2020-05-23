// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test ColorEnriched::internal::find_connection_between_subdomains(...)
// function. Check if the function correctly finds if two subdomains
// share an edge/node.

#include <deal.II/fe/fe_enriched.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <vector>

#include "../tests.h"


/*
 * Construct a class template which finds if a cell
 * is within region or not, based on distance of cell
 * center from region center.
 */
template <int dim>
struct predicate_template
{
  predicate_template(const Point<dim> p, const int radius)
    : p(p)
    , radius(radius)
  {}

  template <class Iterator>
  bool
  operator()(const Iterator &i)
  {
    return ((i->center() - p).norm() < radius);
  }

private:
  Point<dim> p;
  int        radius;
};



template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  // Construct grid
  Triangulation<dim>  triangulation;
  hp::DoFHandler<dim> dof_handler(triangulation);
  GridGenerator::hyper_cube(triangulation, -20, 20);
  triangulation.refine_global(4);

  // Construct vector of predicates for 2 and 3 dimensions
  Assert(dim == 2 || dim == 3, ExcDimensionMismatch2(dim, 2, 3));
  typedef std::function<bool(
    const typename Triangulation<dim>::active_cell_iterator &)>
                                  predicate_function;
  std::vector<predicate_function> predicates;
  predicates.resize(5);
  if (dim == 2)
    {
      // Radius set such that every region has 2^2 = 4 cells
      predicates[1] = predicate_template<dim>(Point<dim>(7.5, 7.5), 2);
      predicates[2] = predicate_template<dim>(Point<dim>(5, 5), 2);
      predicates[0] = predicate_template<dim>(Point<dim>(), 2);
      predicates[3] = predicate_template<dim>(Point<dim>(-5, -5), 2);
      predicates[4] = predicate_template<dim>(Point<dim>(-10, -10), 2);
    }
  else if (dim == 3)
    {
      // Radius set such that every region has 2^3 = 8 cells
      predicates[1] = predicate_template<dim>(Point<dim>(7.5, 7.5, 7.5), 3);
      predicates[2] = predicate_template<dim>(Point<dim>(5, 5, 5), 3);
      predicates[0] = predicate_template<dim>(Point<dim>(), 3);
      predicates[3] = predicate_template<dim>(Point<dim>(-5, -5, -5), 3);
      predicates[4] = predicate_template<dim>(Point<dim>(-10, -10, -10), 3);
    }

  // Check pair-wise connections between predicate regions
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      {
        deallog << i << ":" << j << "="
                << ColorEnriched::internal::find_connection_between_subdomains(
                     dof_handler, predicates[i], predicates[j])
                << std::endl;
      }
}


int
main()
{
  initlog();

  try
    {
      test<2>();
      test<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
  return 0;
}
