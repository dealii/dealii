// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test the rank-2 tensor transformations of MappingFEField by comparing
// against MappingQ. MappingFEField represents the same affine geometry plus
// a translation, so both mappings have the same nontrivial Jacobian.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <string>
#include <vector>

#include "../tests.h"

using namespace dealii;


// Access the mapping data stored inside FEValues.
template <int dim>
class FEValues2 : public FEValues<dim>
{
public:
  using FEValues<dim>::FEValues;

  const typename Mapping<dim>::InternalDataBase &
  get_mapping_data() const
  {
    Assert(this->mapping_data, ExcInternalError());
    return *this->mapping_data;
  }
};


template <int dim>
class Translation : public Function<dim>
{
public:
  Translation()
    : Function<dim>(dim)
  {}

  void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    for (unsigned int d = 0; d < dim; ++d)
      values[d] = p[d] + 0.2 * (d + 1);
  }
};


template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  // Use a nontrivial affine geometry so that the gradient transformations
  // involve a non-identity Jacobian.
  GridTools::transform(
    [](const Point<dim> &p)
    {
      Point<dim> x;

      for (unsigned int d = 0; d < dim; ++d)
        x[d] = (1.2 + 0.3 * d) * p[d];

      if constexpr (dim > 1)
        x[0] += 0.25 * p[1];

      if constexpr (dim > 2)
        x[1] += 0.20 * p[2];

      return x;
    },
    triangulation);

  FE_Q<dim>     fe_q(1);
  FESystem<dim> position_fe(fe_q, dim);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(position_fe);

  Vector<double> position(dof_handler.n_dofs());

  // Add only a translation to the affine geometry. This leaves the Jacobian
  // unchanged, so MappingQ and MappingFEField should transform tensors equally.
  VectorTools::interpolate(dof_handler, Translation<dim>(), position);

  MappingFEField<dim> mapping_fe_field(dof_handler, position);
  MappingQ<dim>       mapping_q(1);

  const QGauss<dim> quadrature(2);

  const UpdateFlags flags =
    update_covariant_transformation | update_contravariant_transformation;

  FEValues2<dim> fe_values_q(mapping_q,
                             position_fe,
                             quadrature,
                             flags);

  FEValues2<dim> fe_values_fe_field(mapping_fe_field,
                                    position_fe,
                                    quadrature,
                                    flags);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values_q.reinit(cell);
      fe_values_fe_field.reinit(cell);

      std::vector<Tensor<2, dim>> input(quadrature.size());

      // Use nonsymmetric input tensors so that index or transpose errors are
      // visible in the comparison.
      for (unsigned int q = 0; q < input.size(); ++q)
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            input[q][i][j] = 1.0 + q + 0.1 * i + 0.01 * j;

      const auto test_mapping =
        [&](const MappingKind mapping_kind, const std::string &name)
        {
          std::vector<Tensor<2, dim>> output_q(quadrature.size());
          std::vector<Tensor<2, dim>> output_fe_field(quadrature.size());

          mapping_q.transform(make_array_view(input),
                              mapping_kind,
                              fe_values_q.get_mapping_data(),
                              make_array_view(output_q));

          mapping_fe_field.transform(
            make_array_view(input),
            mapping_kind,
            fe_values_fe_field.get_mapping_data(),
            make_array_view(output_fe_field));

          for (unsigned int q = 0; q < input.size(); ++q)
            AssertThrow((output_q[q] - output_fe_field[q]).norm() < 1e-12,
                        ExcInternalError());

          deallog << name << ": OK" << std::endl;
        };

      test_mapping(mapping_covariant_gradient,
                   "mapping_covariant_gradient");

      test_mapping(mapping_contravariant_gradient,
                   "mapping_contravariant_gradient");
    }
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
}
