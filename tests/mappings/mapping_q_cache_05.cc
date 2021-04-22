// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// Test MappingQCache initialization with point lambda/Function

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_cache.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

template <int dim>
class Solution : public Function<dim>
{
public:
  Solution(const bool is_displacement_function)
    : Function<dim>(dim)
    , is_displacement_function(is_displacement_function)
  {}

  double
  value(const Point<dim> &point, const unsigned int compontent) const
  {
    return std::sin(point[compontent] * 0.5 * numbers::PI) -
           (is_displacement_function ? point[compontent] : 0.0);
  }

private:
  const bool is_displacement_function;
};

static int counter = 0;

template <int dim, typename Fu>
void
do_test(const unsigned int degree,
        const Fu &         fu,
        const bool         is_displacement_function)
{
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 4);

  MappingQGeneric<dim> mapping(degree);
  MappingQCache<dim>   mapping_cache(degree);
  mapping_cache.initialize(mapping, tria, fu, is_displacement_function);

  {
    DataOut<dim> data_out;

    data_out.attach_triangulation(tria);

    data_out.build_patches(mapping_cache,
                           2,
                           DataOut<dim>::CurvedCellRegion::curved_inner_cells);

#if false
    std::ofstream output("test." + std::to_string(counter++) + ".vtk");
    data_out.write_vtk(output);
#else
    data_out.write_vtk(deallog.get_file_stream());
#endif
  }
}


int
main()
{
  initlog();
  do_test<2>(3, Solution<2>(true), true);
  do_test<2>(3, Solution<2>(false), false);
  do_test<2>(3,
             [](const typename Triangulation<2>::cell_iterator &,
                const Point<2> &p) -> Point<2> {
               Point<2> result;

               for (unsigned int compontent = 0; compontent < 2; ++compontent)
                 result[compontent] =
                   std::sin(p[compontent] * 0.5 * numbers::PI) - p[compontent];

               return result;
             },
             true);
  do_test<2>(3,
             [](const typename Triangulation<2>::cell_iterator &,
                const Point<2> &p) -> Point<2> {
               Point<2> result;

               for (unsigned int compontent = 0; compontent < 2; ++compontent)
                 result[compontent] =
                   std::sin(p[compontent] * 0.5 * numbers::PI);

               return result;
             },
             false);
}
