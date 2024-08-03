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
// FEValuesViews::RenumberedView. Test get_function_values and
// get_function_values_from_local_dof_values

#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_coupling_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

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

  // Create a non trivial function of x and y
  FunctionParser<dim> f("x+10*y");

  Vector<double> vec(dofh.n_dofs());
  VectorTools::interpolate(dofh, f, vec);

  MappingQ<dim> mapping(1);
  UpdateFlags   update_flags =
    update_values | update_quadrature_points | update_JxW_values;

  FEValues<dim> fe_v(mapping, fe, QGauss<dim>(fe.degree + 1), update_flags);
  const auto   &cell = dofh.begin_active();
  fe_v.reinit(cell);

  const auto &fv1_scalar = fe_v[FEValuesExtractors::Scalar(0)];

  // Small print function
  auto print = [](const auto &str, const auto &v) {
    deallog << str << " = ";
    for (const auto &i : v)
      deallog << i << " ";
    deallog << std::endl;
  };

  // Now build an artificial renumbering vector for dofs and one for quadratures
  // - The first renumbering, simply shifts the dof indices by 1, as if there
  //   was an additional dof at the beginning, which is ignored by this view
  // - The second renumbering, repeats each quadrature point twice, so that the
  //   returned values are actually repeated.
  using iota = std_cxx20::ranges::iota_view<unsigned int, unsigned int>;

  std::vector<unsigned int> dof_renumbering(fe_v.dofs_per_cell + 1);
  std::vector<unsigned int> quad_renumbering(fe_v.n_quadrature_points * 2);
  {
    auto id            = iota(0, fe_v.dofs_per_cell);
    dof_renumbering[0] = numbers::invalid_unsigned_int;
    std::copy(id.begin(), id.end(), dof_renumbering.begin() + 1);
  }
  {
    auto id = iota(0, fe_v.n_quadrature_points);
    std::copy(id.begin(), id.end(), quad_renumbering.begin());
    std::copy(id.begin(), id.end(), quad_renumbering.begin() + id.size());
  }

  deallog << "dof_renumbering = ";
  for (const auto &i : dof_renumbering)
    deallog << (int)i << " ";
  deallog << std::endl;
  print("quad_renumbering", quad_renumbering);

  std::vector<double> function_values(fe_v.n_quadrature_points);
  std::vector<double> function_values_renumbered(fe_v.n_quadrature_points * 2);

  fv1_scalar.get_function_values(vec, function_values);
  print("function_values", function_values);

  FEValuesViews::RenumberingData data(fe_v.dofs_per_cell,
                                      fe_v.n_quadrature_points,
                                      dof_renumbering,
                                      quad_renumbering);
  FEValuesViews::RenumberedView<FEValuesViews::Scalar<dim>>
    fv1_scalar_renumbered(fv1_scalar, data);

  fv1_scalar_renumbered.get_function_values(vec, function_values_renumbered);
  print("function_values_renumbered", function_values_renumbered);

  std::vector<double> dof_values(fe_v.dofs_per_cell);
  std::vector<double> dof_values_renumbered(fe_v.dofs_per_cell + 1);

  cell->get_dof_values(vec, dof_values.begin(), dof_values.end());

  // Now fill the first dof with a dummy value
  dof_values_renumbered[0] = 42;
  std::copy(dof_values.begin(),
            dof_values.end(),
            dof_values_renumbered.begin() + 1);

  print("dof_values", dof_values);
  print("dof_values_renumbered", dof_values_renumbered);

  fv1_scalar.get_function_values_from_local_dof_values(dof_values,
                                                       function_values);
  print("function_values from dof_values", function_values);

  fv1_scalar_renumbered.get_function_values_from_local_dof_values(
    dof_values_renumbered, function_values_renumbered);

  print("function_values_renumbered from dof_values_renumbered",
        function_values_renumbered);
}

int
main()
{
  initlog();
  test<2>();
}
