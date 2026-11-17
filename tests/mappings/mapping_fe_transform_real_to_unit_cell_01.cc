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


// Test MappingFE::transform_real_to_unit_cell() on curved quadrilaterals,
// triangles, hexahedra, and tetrahedra: map the quadrature points of the cell
// and of its faces to real space with transform_unit_to_real_cell() and back,
// which has to reproduce the points. The meshes of the unit square and cube
// are mapped to a sector of a cylinder shell with a manifold on all cells,
// i.e., all cells are curved.

#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


// Map the unit square or cube to a sector of a cylinder shell around the
// z-axis with the radii 1 and 2 and the opening angle pi/2: the first
// coordinate is the radius, the second one the angle.
template <int dim>
class ShellSector : public ChartManifold<dim>
{
public:
  std::unique_ptr<Manifold<dim>>
  clone() const override
  {
    return std::make_unique<ShellSector<dim>>();
  }

  Point<dim>
  push_forward(const Point<dim> &chart_point) const override
  {
    Point<dim>   p     = chart_point;
    const double angle = 0.5 * numbers::PI * chart_point[1];
    p[0]               = (1. + chart_point[0]) * std::cos(angle);
    p[1]               = (1. + chart_point[0]) * std::sin(angle);
    return p;
  }

  Point<dim>
  pull_back(const Point<dim> &p) const override
  {
    Point<dim> chart_point = p;
    chart_point[0]         = std::sqrt(p[0] * p[0] + p[1] * p[1]) - 1.;
    chart_point[1]         = std::atan2(p[1], p[0]) / (0.5 * numbers::PI);
    return chart_point;
  }
};



template <int dim>
void
map_to_shell_sector(Triangulation<dim> &tria)
{
  const ShellSector<dim> manifold;
  GridTools::transform(
    [&](const Point<dim> &p) { return manifold.push_forward(p); }, tria);
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);
}



template <int dim>
void
test(const Triangulation<dim> &tria, const FiniteElement<dim> &fe_mapping)
{
  const MappingFE<dim>     mapping(fe_mapping);
  const ReferenceCell<dim> reference_cell = fe_mapping.reference_cell();

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

  // All cells have to be valid, i.e., the map is invertible.
  const FE_Nothing<dim> fe(reference_cell);
  FEValues<dim>         fe_values(mapping,
                          fe,
                          reference_cell.get_gauss_type_quadrature(5),
                          update_jacobians);

  unsigned int n_points = 0, n_failed = 0;
  for (const auto &cell : tria.active_cell_iterators())
    {
      fe_values.reinit(cell);
      for (const unsigned int q : fe_values.quadrature_point_indices())
        AssertThrow(fe_values.jacobian(q).determinant() > 0,
                    ExcInternalError());

      for (const Point<dim> &p_unit : unit_points)
        {
          ++n_points;
          const Point<dim> p_real =
            mapping.transform_unit_to_real_cell(cell, p_unit);
          try
            {
              n_failed += mapping.transform_real_to_unit_cell(cell, p_real)
                            .distance(p_unit) > 1e-8;
            }
          catch (const typename Mapping<dim>::ExcTransformationFailed &)
            {
              ++n_failed;
            }
        }
    }

  deallog << fe_mapping.get_name() << ": failed round trips: " << n_failed
          << " of " << n_points << " points" << std::endl;
}



int
main()
{
  initlog();

  {
    Triangulation<2> tria;
    GridGenerator::subdivided_hyper_cube(tria, 2);
    map_to_shell_sector(tria);
    for (unsigned int degree = 1; degree <= 3; ++degree)
      test(tria, FE_Q<2>(degree));
  }
  {
    Triangulation<2> tria;
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 2);
    map_to_shell_sector(tria);
    for (unsigned int degree = 1; degree <= 3; ++degree)
      test(tria, FE_SimplexP<2>(degree));
  }
  {
    Triangulation<3> tria;
    GridGenerator::subdivided_hyper_cube(tria, 2);
    map_to_shell_sector(tria);
    for (unsigned int degree = 1; degree <= 3; ++degree)
      test(tria, FE_Q<3>(degree));
  }
  {
    Triangulation<3> tria;
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 2);
    map_to_shell_sector(tria);
    for (unsigned int degree = 1; degree <= 3; ++degree)
      test(tria, FE_SimplexP<3>(degree));
  }
}
