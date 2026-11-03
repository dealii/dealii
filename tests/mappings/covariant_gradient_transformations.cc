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

// Check covariant gradient transformations against C G C^T for a known affine
// geometry, including nonsymmetric tensors and rectangular Jacobians.

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_p1.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools_interpolate.h>

#include <string>
#include <vector>

#include "../tests.h"


// Access the mapping data initialized by FEValues::reinit().
template <int dim, int spacedim>
class TestFEValues : public FEValues<dim, spacedim>
{
public:
  using FEValues<dim, spacedim>::FEValues;

  const typename Mapping<dim, spacedim>::InternalDataBase &
  get_mapping_data() const
  {
    return *this->mapping_data;
  }
};


template <int dim, int spacedim>
void
check(const Mapping<dim, spacedim>           &mapping,
      const FiniteElement<dim, spacedim>     &fe,
      const Quadrature<dim>                  &quadrature,
      const Triangulation<dim, spacedim>     &tria,
      const DerivativeForm<1, dim, spacedim> &covariant,
      const std::string                      &name)
{
  TestFEValues<dim, spacedim> fe_values(mapping,
                                        fe,
                                        quadrature,
                                        update_covariant_transformation |
                                          update_contravariant_transformation);

  std::vector<Tensor<2, dim>>      input(quadrature.size());
  std::vector<Tensor<2, spacedim>> output(quadrature.size());
  for (unsigned int q = 0; q < input.size(); ++q)
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        input[q][i][j] = 1.0 + 0.3 * q + 0.2 * i - 0.1 * j;

  for (const auto &cell : tria.active_cell_iterators())
    {
      fe_values.reinit(cell);
      mapping.transform(make_array_view(input),
                        mapping_covariant_gradient,
                        fe_values.get_mapping_data(),
                        make_array_view(output));

      for (unsigned int q = 0; q < input.size(); ++q)
        {
          Tensor<2, spacedim> expected;
          for (unsigned int i = 0; i < spacedim; ++i)
            for (unsigned int j = 0; j < spacedim; ++j)
              for (unsigned int k = 0; k < dim; ++k)
                for (unsigned int l = 0; l < dim; ++l)
                  expected[i][j] +=
                    covariant[i][k] * input[q][k][l] * covariant[j][l];
          AssertThrow((output[q] - expected).norm() < 1e-11,
                      ExcMessage(name + ": incorrect covariant gradient"));
        }
    }
  deallog << name << ": OK" << std::endl;
}


template <int dim, int spacedim>
void
test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  // This fixed Jacobian has unequal scales, skew, and nonzero embedding
  // components. Compute C = J (J^T J)^-1 independently of the mapping.
  DerivativeForm<1, dim, spacedim> jacobian;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      jacobian[i][j] =
        i == j ? 1.5 + 0.25 * i : 0.07 * (i + 1) + 0.11 * (j + 1);

  Tensor<2, dim> gram;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < spacedim; ++k)
        gram[i][j] += jacobian[k][i] * jacobian[k][j];
  const Tensor<2, dim> inverse_gram = invert(gram);

  DerivativeForm<1, dim, spacedim> covariant;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        covariant[i][j] += jacobian[i][k] * inverse_gram[k][j];

  const auto transform = [&jacobian](const Point<spacedim> &p) {
    Point<spacedim> result;
    for (unsigned int i = 0; i < spacedim; ++i)
      {
        result[i] = 0.2 * (i + 1);
        for (unsigned int j = 0; j < dim; ++j)
          result[i] += jacobian[i][j] * p[j];
      }
    return result;
  };

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  GridTools::transform(transform, tria);

  FE_Q<dim, spacedim>            fe(1);
  QGauss<dim>                    quadrature(2);
  MappingQ<dim, spacedim>        mapping_q(1);
  MappingManifold<dim, spacedim> mapping_manifold;
  check(mapping_q, fe, quadrature, tria, covariant, "MappingQ");
  check(mapping_manifold, fe, quadrature, tria, covariant, "MappingManifold");

  FESystem<dim, spacedim>   position_fe(fe, spacedim);
  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(position_fe);
  Vector<double> position(dof_handler.n_dofs());
  VectorTools::get_position_vector(dof_handler, position);
  MappingFEField<dim, spacedim> mapping_fe_field(dof_handler, position);
  check(mapping_fe_field, fe, quadrature, tria, covariant, "MappingFEField");

  Triangulation<dim, spacedim> simplex_tria;
  GridGenerator::reference_cell(simplex_tria,
                                ReferenceCells::get_simplex<dim>());
  GridTools::transform(transform, simplex_tria);
  FE_SimplexP<dim, spacedim> simplex_fe(1);
  QGaussSimplex<dim>         simplex_quadrature(2);
  MappingFE<dim, spacedim>   mapping_fe(simplex_fe);
  MappingP1<dim, spacedim>   mapping_p1;
  check(mapping_fe,
        simplex_fe,
        simplex_quadrature,
        simplex_tria,
        covariant,
        "MappingFE");
  check(mapping_p1,
        simplex_fe,
        simplex_quadrature,
        simplex_tria,
        covariant,
        "MappingP1");
}


int
main()
{
  initlog();
  test<1, 1>();
  test<2, 2>();
  test<3, 3>();
  test<1, 2>();
  test<1, 3>();
  test<2, 3>();
}
