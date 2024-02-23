// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test created to reproduce a bug in make_flux_sparsity_pattern that occurs
// in 1D when neighboring cells are on different levels. Set up a triangulation
// with 3 cells in 1D, where the left-most element is on level 0 and the middle
// and right-most element is on level 1:
//
// |--l0--|-l1-|-l1-|.
//
// Distribute FE_Q<1>(1) elements, call make_flux_sparsity_pattern and check
// that the created pattern is the expected one.


#include <deal.II/base/table.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


/*
 * Call the version of make_flux_sparsity_pattern that takes both cell and face
 * integral couplings. Print the constructed sparsity pattern to deallog.
 */
template <int dim>
void
create_and_print_flux_sparsity_pattern(const DoFHandler<dim> &dof_handler)
{
  AffineConstraints<double> constraints;
  constraints.close();

  const bool keep_constrained_dofs = true;

  // Set all couplings to be connected.
  const unsigned int           n_components = 1;
  Table<2, DoFTools::Coupling> couplings(n_components, n_components);
  couplings[0][0] = DoFTools::Coupling::always;
  Table<2, DoFTools::Coupling> face_couplings(n_components, n_components);
  face_couplings[0][0] = DoFTools::Coupling::always;

  const types::subdomain_id subdomain_id = 0;

  DynamicSparsityPattern dynamic_pattern(dof_handler.n_dofs());

  DoFTools::make_flux_sparsity_pattern(dof_handler,
                                       dynamic_pattern,
                                       constraints,
                                       keep_constrained_dofs,
                                       couplings,
                                       face_couplings,
                                       subdomain_id);

  dynamic_pattern.print(deallog.get_file_stream());
}



/*
 * Create the 3-cells-triangulation (described at the top) by first creating 2
 * cells one the same level and then refining the right one.
 */
void
create_3_elements_on_2_different_levels(Triangulation<1> &triangulation)
{
  const unsigned int n_elements = 3;
  GridGenerator::subdivided_hyper_cube(triangulation, n_elements - 1);
  auto cell = triangulation.begin_active();
  cell++;
  cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
}



int
main()
{
  initlog();
  const int dim = 1;

  Triangulation<dim> triangulation;
  create_3_elements_on_2_different_levels(triangulation);

  const FE_Q<dim> element(1);

  // Since the implementation of make_flux_sparsity_pattern is specialized for
  // hp::DoFHandler we create and print the sparsity pattern for both
  // DoFHandler-types.
  deallog << "dealii::DoFHandler" << std::endl;
  {
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(element);
    create_and_print_flux_sparsity_pattern(dof_handler);
  }
  deallog << "hp::DoFHandler" << std::endl;
  {
    hp::FECollection<dim> fe_collection(element);
    DoFHandler<dim>       dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_collection);
    create_and_print_flux_sparsity_pattern(dof_handler);
  }
}
