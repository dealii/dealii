// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check FESystem system_to_component_index() and system_to_base_index()
// for a simple FE_Q case and prints the results.
// If at some point these results will change (due to internal refactoring),
// update the documentation of FiniteElement::system_to_base_index()

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <iomanip>
#include <string>

#include "../tests.h"

template <int dim>
void
test(const bool renumber = false)
{
  Triangulation<dim> triangulation;
  FESystem<dim>      fe_basis(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
  DoFHandler<dim>    dof_handler(triangulation);
  GridGenerator::hyper_cube(triangulation);
  dof_handler.distribute_dofs(fe_basis);

  if (renumber)
    {
      std::vector<unsigned int> component_to_block_indices(dim + 1);
      for (int i = 0; i < dim; ++i)
        component_to_block_indices[i] = 0;
      component_to_block_indices[dim] = 1;
      DoFRenumbering::component_wise(dof_handler, component_to_block_indices);
    }

  const unsigned int n_dofs = dof_handler.n_dofs();

  std::cout
    << " * | DoF    | Component  | Base element | Shape function within base | Multiplicity |"
    << std::endl
    << " * | :----: | :--------: | :----------: | :------------------------: | :----------: |"
    << std::endl;

  for (unsigned int i = 0; i < n_dofs; ++i)
    {
      const unsigned int component =
        fe_basis.system_to_component_index(i).first;
      const unsigned int within_base =
        fe_basis.system_to_component_index(i).second;
      const unsigned int base = fe_basis.system_to_base_index(i).first.first;
      const unsigned int multiplicity =
        fe_basis.system_to_base_index(i).first.second;
      const unsigned int within_base_ =
        fe_basis.system_to_base_index(i).second; // same as above
      std::cout << std::setfill(' ') << " * | " << std::setw(6) << i << " | "
                << std::setw(10) << component << " | " << std::setw(12) << base
                << " | " << std::setw(26) << within_base << " | "
                << std::setw(12) << multiplicity << " |" << std::endl;
    }

  // print grid and DoFs for visual inspection
  if (true)
    {
      std::map<types::global_dof_index, Point<dim>> support_points;
      MappingQ1<dim>                                mapping;
      DoFTools::map_dofs_to_support_points(mapping,
                                           dof_handler,
                                           support_points);

      const std::string filename = "grid" + Utilities::int_to_string(dim) +
                                   Utilities::int_to_string(renumber) + ".gp";
      std::ofstream f(filename);

      f << "set terminal png size 420,440 enhanced font \"Helvetica,16\""
        << std::endl
        << "set output \"grid" << Utilities::int_to_string(dim)
        << Utilities::int_to_string(renumber) << ".png\"" << std::endl
        << "set size square" << std::endl
        << "set view equal xy" << std::endl
        << "unset xtics" << std::endl
        << "unset ytics" << std::endl
        << "unset border" << std::endl
        << "set xrange [0: 1.05]" << std::endl
        << "set yrange [0: 1.05]" << std::endl
        << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 0.5,0.5 notitle"
        << std::endl;
      GridOut().write_gnuplot(triangulation, f);
      f << 'e' << std::endl;

      DoFTools::write_gnuplot_dof_support_point_info(f, support_points);
      f << 'e' << std::endl;
    }
}

int
main()
{
  test<2>(false);

  return 0;
}
