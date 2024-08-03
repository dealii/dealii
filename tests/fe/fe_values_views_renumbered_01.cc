// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2024 by the deal.II authors
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


// Test renumbering of FEValuesViews::Scalar using
// FEValuesViews::RenumberedView

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_coupling_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

#include "../test_grids.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  DoFHandler<dim> dofh(tria);
  FE_Q<dim>       fe(1);
  dofh.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  UpdateFlags   update_flags =
    update_values | update_quadrature_points | update_JxW_values;

  FEValues<dim> fv1(mapping, fe, QGauss<dim>(fe.degree + 1), update_flags);
  fv1.reinit(dofh.begin_active());

  FEValuesViews::Scalar<dim> fv1_scalar(fv1, 0);

  using iota = std_cxx20::ranges::iota_view<unsigned int, unsigned int>;
  {
    // No renumbering
    for (const auto &i : fv1.dof_indices())
      {
        deallog << "i = " << i << ", fv1_scalar = " << fv1_scalar.value(i, 0)
                << std::endl;
      }
  }
  {
    // Trivial renumbering
    auto                           id = iota(0, fv1.dofs_per_cell);
    std::vector<unsigned int>      renumber(id.begin(), id.end());
    FEValuesViews::RenumberingData data(fv1.dofs_per_cell,
                                        fv1.n_quadrature_points,
                                        renumber);
    FEValuesViews::RenumberedView<FEValuesViews::Scalar<dim>>
      fv1_scalar_reordered(fv1_scalar, data);

    for (const auto &i : fv1.dof_indices())
      {
        deallog << "i = " << i << ", fv1_scalar_trivial_renumber = "
                << fv1_scalar_reordered.value(i, 0) << std::endl;
      }
  }
  {
    // Inverse renumbering
    auto                      id = iota(0, fv1.dofs_per_cell);
    std::vector<unsigned int> renumber(id.begin(), id.end());
    std::reverse(renumber.begin(), renumber.end());
    FEValuesViews::RenumberingData data(fv1.dofs_per_cell,
                                        fv1.n_quadrature_points,
                                        renumber);
    FEValuesViews::RenumberedView<FEValuesViews::Scalar<dim>>
      fv1_scalar_reordered(fv1_scalar, data);
    for (const auto &i : fv1.dof_indices())
      {
        deallog << "i = " << i << ", fv1_scalar_inverse_renumber = "
                << fv1_scalar_reordered.value(i, 0) << std::endl;
      }
  }
  // Now test quadrature renumbering
  {
    // No renumbering
    for (const auto &q : fv1.quadrature_point_indices())
      {
        deallog << "q = " << q
                << ", fv1_scalar_0_q = " << fv1_scalar.value(0, q) << std::endl;
      }
  }
  {
    // Trivial renumbering
    auto                           id = iota(0, fv1.n_quadrature_points);
    std::vector<unsigned int>      renumber(id.begin(), id.end());
    FEValuesViews::RenumberingData data(fv1.dofs_per_cell,
                                        fv1.n_quadrature_points,
                                        {},
                                        renumber);
    FEValuesViews::RenumberedView<FEValuesViews::Scalar<dim>>
      fv1_scalar_reordered(fv1_scalar, data);
    for (const auto &q : fv1.quadrature_point_indices())
      {
        deallog << "q = " << q << ", fv1_scalar_0_q_trivial_renumbering = "
                << fv1_scalar_reordered.value(0, q) << std::endl;
      }
  }
  {
    // Inverse renumbering
    auto                      id = iota(0, fv1.n_quadrature_points);
    std::vector<unsigned int> renumber(id.begin(), id.end());
    std::reverse(renumber.begin(), renumber.end());
    FEValuesViews::RenumberingData data(fv1.dofs_per_cell,
                                        fv1.n_quadrature_points,
                                        {},
                                        renumber);
    FEValuesViews::RenumberedView<FEValuesViews::Scalar<dim>>
      fv1_scalar_reordered(fv1_scalar, data);
    for (const auto &q : fv1.quadrature_point_indices())
      {
        deallog << "q = " << q << ", fv1_scalar_0_q_inverse_renumbering = "
                << fv1_scalar_reordered.value(0, q) << std::endl;
      }
  }
}

int
main()
{
  initlog();
  test<2>();
}
