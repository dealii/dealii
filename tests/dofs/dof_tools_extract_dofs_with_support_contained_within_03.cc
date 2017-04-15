// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
//

// test DoFTools::extract_dofs_with_support_contained_within() but
// for a slightly different configuration

#include "../tests.h"
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include <set>
#include <cstdio>


template<int dim>
bool
pred_left(const typename Triangulation<dim>::active_cell_iterator &cell)
{
  return (cell->center()(0) < 0.49);
}

template<int dim>
bool
pred_right(const typename Triangulation<dim>::active_cell_iterator &cell)
{
  return (cell->center()(0) > 0.51);
}


template<int dim>
bool
pred_r(const typename Triangulation<dim>::active_cell_iterator &cell)
{
  return (cell->center()(0) < 0.49 &&
          cell->center()(1) < 0.49) ||
         (cell->center()(0) > 0.49 &&
          cell->center()(1) > 0.49);
}



template <int dim>
void test (const bool left = true)
{

  // Setup system
  Triangulation<dim> triangulation;

  GridGenerator::hyper_rectangle (triangulation,
                                  Point<dim>(0,0),
                                  Point<dim>(1,1));

  triangulation.refine_global(1);

  // Extra refinement to generate hanging nodes
  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (pred_r<dim>(cell))
      cell->set_refine_flag ();

  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement ();

  DoFHandler<dim> dh (triangulation);

  FE_Q<dim> fe(2);
  dh.distribute_dofs (fe);

  ConstraintMatrix cm;
  DoFTools::make_hanging_node_constraints(dh,cm);

  IndexSet support = left ?
                     DoFTools::extract_dofs_with_support_contained_within(dh, pred_left<dim>, cm) :
                     DoFTools::extract_dofs_with_support_contained_within(dh, pred_right<dim>, cm);
  support.print(deallog);

  // print grid and DoFs for visual inspection
  if (false)
    {
      std::cout << "-------------------- " << Utilities::int_to_string(dim) << std::endl;
      cm.print(std::cout);

      std::map<types::global_dof_index, Point<dim> > support_points;
      MappingQ1<dim> mapping;
      DoFTools::map_dofs_to_support_points(mapping, dh, support_points);

      const std::string filename =
        "grid" + Utilities::int_to_string(dim) + ".gp";
      std::ofstream f(filename.c_str());

      f << "set terminal png size 400,410 enhanced font \"Helvetica,8\"" << std::endl
        << "set output \"grid" << Utilities::int_to_string(dim) << ".png\"" << std::endl
        << "set size square" << std::endl
        << "set view equal xy" << std::endl
        << "unset xtics" << std::endl
        << "unset ytics" << std::endl
        << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 1,1 notitle" << std::endl;
      GridOut().write_gnuplot (triangulation, f);
      f << "e" << std::endl;

      DoFTools::write_gnuplot_dof_support_point_info(f,
                                                     support_points);
      f << "e" << std::endl;
    }

  dh.clear();
}

int main(int argc, char **argv)
{
  initlog();

  test<2>(true);
  test<2>(false);

  return 0;
}
