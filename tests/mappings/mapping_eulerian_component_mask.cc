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


// Select displacement components from a mixed field, including in codimension.
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1_eulerian.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// Every component has a different polynomial. Thus selecting the wrong field
// cannot produce the expected displacement or its derivative.
template <int spacedim>
class Field : public Function<spacedim>
{
public:
  Field(const unsigned int n_components, const unsigned int degree)
    : Function<spacedim>(n_components)
    , degree(degree)
  {}

  double
  value(const Point<spacedim> &p, const unsigned int c = 0) const override
  {
    return 0.02 * (c + 1) + 0.01 * (c + 2) * p[0] - 0.005 * (c + 1) * p[1] +
           (degree == 2 ? 0.003 * (c + 1) * p[0] * p[0] : 0.);
  }

  double
  derivative(const Point<spacedim> &p,
             const unsigned int     c,
             const unsigned int     d) const
  {
    if (d == 0)
      return 0.01 * (c + 2) + (degree == 2 ? 0.006 * (c + 1) * p[0] : 0.);
    return d == 1 ? -0.005 * (c + 1) : 0.;
  }

private:
  const unsigned int degree;
};

template <int dim, int spacedim>
void
check(const Mapping<dim, spacedim>             &mapping,
      const Triangulation<dim, spacedim>       &tria,
      const Field<spacedim>                    &field,
      const std::array<unsigned int, spacedim> &components)
{
  std::vector<Point<dim>> points;
  for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
    points.push_back(GeometryInfo<dim>::unit_cell_vertex(v));
  for (unsigned int q = 0; q < 3; ++q)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = 0.17 + 0.11 * d + 0.13 * q;
      points.push_back(p);
    }
  const Quadrature<dim>   quadrature(points);
  FE_Q<dim, spacedim>     fe(1);
  FEValues<dim, spacedim> values(mapping, fe, quadrature, update_jacobians);
  values.reinit(tria.begin_active());
  for (unsigned int q = 0; q < points.size(); ++q)
    {
      Point<spacedim> reference;
      for (unsigned int d = 0; d < dim; ++d)
        reference[d] = points[q][d];
      Point<spacedim> expected = reference;
      for (unsigned int c = 0; c < spacedim; ++c)
        expected[c] += field.value(reference, components[c]);
      AssertThrow(
        mapping.transform_unit_to_real_cell(tria.begin_active(), points[q])
            .distance(expected) < 1e-12,
        ExcInternalError());
      for (unsigned int c = 0; c < spacedim; ++c)
        for (unsigned int d = 0; d < dim; ++d)
          AssertThrow(
            std::abs(values.jacobian(q)[c][d] -
                     ((c == d ? 1. : 0.) +
                      field.derivative(reference, components[c], d))) < 1e-12,
            ExcInternalError());
    }
}

template <typename Callable>
void
must_throw(const Callable &call)
{
  bool caught = false;
  try
    {
      call();
    }
  catch (const ExceptionBase &)
    {
      caught = true;
    }
  AssertThrow(caught, ExcInternalError());
}

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  std::array<unsigned int, spacedim> selected, first;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      selected[d] = 2 * d + 1;
      first[d]    = d;
    }
  for (unsigned int degree = 1; degree <= 2; ++degree)
    {
      // The unrelated DG field has no vertex DoFs and a different degree.
      FESystem<dim, spacedim>   fe(FE_Q<dim, spacedim>(degree),
                                 2 * spacedim + 1,
                                 FE_DGQ<dim, spacedim>(degree + 1),
                                 1);
      DoFHandler<dim, spacedim> dofs(tria);
      Vector<double>            displacement;
      std::vector<bool>         selection(fe.n_components(), false);
      for (auto c : selected)
        selection[c] = true;
      const ComponentMask mask(selection);
      // MappingQ1Eulerian must remain constructible before distribute_dofs().
      MappingQ1Eulerian<dim, Vector<double>, spacedim> q1(dofs,
                                                          displacement,
                                                          mask);
      dofs.distribute_dofs(fe);
      displacement.reinit(dofs.n_dofs());
      Field<spacedim> field(fe.n_components(), degree);
      VectorTools::interpolate(dofs, field, displacement);
      MappingQEulerian<dim, Vector<double>, spacedim> q(degree,
                                                        dofs,
                                                        displacement,
                                                        mask);
      check<dim, spacedim>(q, tria, field, selected);
      check<dim, spacedim>(*q.clone(), tria, field, selected);
      if (degree == 1)
        {
          check<dim, spacedim>(q1, tria, field, selected);
          check<dim, spacedim>(*q1.clone(), tria, field, selected);
        }
      else
        {
          // Higher-order FE_Q is also valid: its nonvertex shapes vanish at
          // vertices. An affine field is reproduced by the Q1 geometry.
          Field<spacedim> affine(fe.n_components(), 1);
          VectorTools::interpolate(dofs, affine, displacement);
          check<dim, spacedim>(q1, tria, affine, selected);
          check<dim, spacedim>(*q1.clone(), tria, affine, selected);
          VectorTools::interpolate(dofs, field, displacement);
        }
      // Both old signatures and an explicitly empty mask use the first fields.
      MappingQEulerian<dim, Vector<double>, spacedim> legacy(
        degree, dofs, displacement, numbers::invalid_unsigned_int);
      MappingQEulerian<dim, Vector<double>, spacedim> empty(degree,
                                                            dofs,
                                                            displacement,
                                                            ComponentMask());
      check<dim, spacedim>(legacy, tria, field, first);
      check<dim, spacedim>(empty, tria, field, first);
    }
  {
    // This element's cell-constant mode contributes at every vertex.
    // Silently ignoring it would produce incorrect displaced geometry.
    FESystem<dim, spacedim>   unsupported(FE_Q_DG0<dim, spacedim>(1), spacedim);
    DoFHandler<dim, spacedim> dofs(tria);
    dofs.distribute_dofs(unsupported);
    Vector<double> displacement(dofs.n_dofs());
    displacement = 1.;
    must_throw([&]() {
      MappingQ1Eulerian<dim, Vector<double>, spacedim> bad(
        dofs, displacement, ComponentMask(std::vector<bool>(spacedim, true)));
      bad.get_vertices(tria.begin_active());
    });
  }
  FESystem<dim, spacedim>   fe(FE_Q<dim, spacedim>(1), spacedim);
  DoFHandler<dim, spacedim> dofs(tria);
  Vector<double>            displacement;
  MappingQ1Eulerian<dim, Vector<double>, spacedim> legacy(dofs, displacement);
  dofs.distribute_dofs(fe);
  displacement.reinit(dofs.n_dofs());
  Field<spacedim> field(spacedim, 1);
  VectorTools::interpolate(dofs, field, displacement);
  check<dim, spacedim>(legacy, tria, field, first);
  deallog << dim << ',' << spacedim << " OK" << std::endl;
}

int
main()
{
  initlog();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
