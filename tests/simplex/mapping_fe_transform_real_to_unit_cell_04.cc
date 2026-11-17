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


// Test MappingFE::transform_real_to_unit_cell() on strongly curved simplices:
// map the quadrature points of the cell and of its faces and the points of a
// lattice of the reference cell to real space with
// transform_unit_to_real_cell() and back, which has to reproduce the points.
// The reference cell is displaced by a sine wave, i.e., the cell is valid but
// far from affine. For about half of the points, the Newton iteration does not
// converge from the inverse of the affine map through the vertices, but from
// the barycenter of the reference cell, which is the fallback. For some points,
// the iteration from the affine map converges to another preimage outside of
// the reference cell, which also requires the fallback.

#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


// Flat manifold displaced by the sine wave
// d_i(x) = amplitude_i * sin(wave_vector_i * x + phase_i).
template <int dim>
class SineDisplacement : public Manifold<dim>
{
public:
  SineDisplacement(const Tensor<1, dim> &amplitude,
                   const Tensor<2, dim> &wave_vector,
                   const Tensor<1, dim> &phase)
    : amplitude(amplitude)
    , wave_vector(wave_vector)
    , phase(phase)
  {}

  std::unique_ptr<Manifold<dim>>
  clone() const override
  {
    return std::make_unique<SineDisplacement<dim>>(*this);
  }

  Point<dim>
  get_new_point(const ArrayView<const Point<dim>> &surrounding_points,
                const ArrayView<const double>     &weights) const override
  {
    Point<dim> p;
    for (unsigned int i = 0; i < surrounding_points.size(); ++i)
      p += weights[i] * surrounding_points[i];

    Point<dim> p_displaced = p;
    for (unsigned int d = 0; d < dim; ++d)
      p_displaced[d] += amplitude[d] * std::sin(wave_vector[d] * p + phase[d]);
    return p_displaced;
  }

private:
  const Tensor<1, dim> amplitude;
  const Tensor<2, dim> wave_vector;
  const Tensor<1, dim> phase;
};



// Points of the reference simplex with the barycentric coordinates i / n,
// i.e., including the vertices, edges, and faces.
template <int dim>
std::vector<Point<dim>>
lattice_points(const unsigned int n)
{
  std::vector<Point<dim>> points;
  if constexpr (dim == 2)
    for (unsigned int j = 0; j <= n; ++j)
      for (unsigned int i = 0; i + j <= n; ++i)
        points.emplace_back(1. * i / n, 1. * j / n);
  else
    for (unsigned int k = 0; k <= n; ++k)
      for (unsigned int j = 0; j + k <= n; ++j)
        for (unsigned int i = 0; i + j + k <= n; ++i)
          points.emplace_back(1. * i / n, 1. * j / n, 1. * k / n);
  return points;
}



template <int dim>
void
test(const Manifold<dim> &manifold, const unsigned int n_lattice)
{
  const ReferenceCell<dim> reference_cell = ReferenceCells::get_simplex<dim>();

  Triangulation<dim> tria;
  GridGenerator::reference_cell(tria, reference_cell);
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);

  const FE_SimplexP<dim> fe_mapping(2);
  const MappingFE<dim>   mapping(fe_mapping);

  std::vector<Point<dim>> unit_points =
    reference_cell.get_gauss_type_quadrature(3).get_points();
  for (const unsigned int f : reference_cell.face_indices())
    {
      const Quadrature<dim> quadrature_on_face =
        QProjector<dim>::project_to_face(
          reference_cell,
          reference_cell.face_reference_cell(f).get_gauss_type_quadrature(3),
          f,
          numbers::default_geometric_orientation);
      unit_points.insert(unit_points.end(),
                         quadrature_on_face.get_points().begin(),
                         quadrature_on_face.get_points().end());
    }

  const std::vector<Point<dim>> lattice = lattice_points<dim>(n_lattice);
  unit_points.insert(unit_points.end(), lattice.begin(), lattice.end());

  // The cell has to be valid, i.e., the map is invertible. The Jacobian is
  // checked at the points of the lattice, which includes the boundary, where
  // the cell is inverted first for larger amplitudes.
  const FE_Nothing<dim> fe(reference_cell);
  FEValues<dim>         fe_values(mapping,
                          fe,
                          Quadrature<dim>(lattice),
                          update_jacobians);

  const auto cell = tria.begin_active();
  fe_values.reinit(cell);
  for (const unsigned int q : fe_values.quadrature_point_indices())
    AssertThrow(fe_values.jacobian(q).determinant() > 0, ExcInternalError());

  unsigned int n_failed = 0;
  for (const Point<dim> &p_unit : unit_points)
    {
      const Point<dim> p_real =
        mapping.transform_unit_to_real_cell(cell, p_unit);
      try
        {
          n_failed +=
            mapping.transform_real_to_unit_cell(cell, p_real).distance(p_unit) >
            1e-8;
        }
      catch (const typename Mapping<dim>::ExcTransformationFailed &)
        {
          ++n_failed;
        }
    }

  deallog << dim << "D " << fe_mapping.get_name()
          << ": failed round trips: " << n_failed << " of "
          << unit_points.size() << " points" << std::endl;
}



int
main()
{
  initlog();

  test(SineDisplacement<2>(Tensor<1, 2>({-0.18, -0.72}),
                           Tensor<2, 2>({{5., -4.}, {-6., 6.}}),
                           Tensor<1, 2>({0., -2.})),
       40);
  test(SineDisplacement<3>(Tensor<1, 3>({0.48, -0.12, 0.48}),
                           Tensor<2, 3>(
                             {{0., -6., 6.}, {-1., -2., -2.}, {1., 1., -5.}}),
                           Tensor<1, 3>({-1., -1., -3.})),
       16);
}
