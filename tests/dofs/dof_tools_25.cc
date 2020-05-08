// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II Authors
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

// Test the version of make_flux_sparsity_pattern that takes an additional
// std::function describing over which faces flux couplings occur.
//
// Set up a 1D-triangulation with 3 cells (left, middle and right) where the
// face between left and middle is located at x=0:
//
// |-L-|-M-|-R-|
//    x=0
//
// Distribute FE_Q<1>(1) elements over this grid. Then call
// make_flux_sparsity_pattern with an std::function specifying that there should
// only be a flux coupling between the left and middle cell. Verify that this
// is the case for the created sparsity pattern.

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


// Call the make_flux_sparsity_pattern with a lambda specifying that only the
// face at x = 0 should have a flux coupling. Print the constructed sparsity
// pattern to deallog.
template <int dim, class DoFHandlerType>
void
create_and_print_pattern(const DoFHandlerType &dof_handler)
{
  DynamicSparsityPattern dynamic_pattern(dof_handler.n_dofs());

  auto face_has_flux_coupling =
    [](const typename DoFHandlerType::active_cell_iterator &cell,
       const unsigned int                                   face_index) {
      // Only add a flux coupling if the face is at x = 0.
      const Point<dim> &center = cell->face(face_index)->center();
      return std::abs(center[0]) < 1e-3;
      return true;
    };

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

  DoFTools::make_flux_sparsity_pattern(dof_handler,
                                       dynamic_pattern,
                                       constraints,
                                       keep_constrained_dofs,
                                       couplings,
                                       face_couplings,
                                       subdomain_id,
                                       face_has_flux_coupling);


  dynamic_pattern.print(deallog.get_file_stream());
}



// The implementation is different depending on which DoFHandler is used. Test
// both of them.
template <int dim>
void
test_with_both_dof_handlers(const Triangulation<dim> &triangulation)
{
  const FE_Q<dim> element(1);

  deallog << "dealii::DoFHandler" << std::endl;
  {
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(element);
    create_and_print_pattern<dim>(dof_handler);
  }

  deallog << std::endl;

  deallog << "hp::DoFHandler" << std::endl;
  {
    hp::FECollection<dim> fe_collection(element);
    hp::DoFHandler<dim>   dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_collection);
    create_and_print_pattern<dim>(dof_handler);
  }
}



//  Create the 3-cell-triangulation described at the top.
void create_3_cell_triangulation(Triangulation<1> &triangulation)
{
  const unsigned int n_elements = 3;
  const double       left       = -1;
  const double       right      = 2;

  GridGenerator::subdivided_hyper_cube(triangulation, n_elements, left, right);
}



int
main()
{
  initlog();

  const int          dim = 1;
  Triangulation<dim> triangulation;
  create_3_cell_triangulation(triangulation);

  test_with_both_dof_handlers(triangulation);
  return 0;
}
