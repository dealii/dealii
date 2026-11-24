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

// Interpolate vector fields with H(div) and H(curl) conforming elements on
// meshes where FEValues changes the sign of some shape functions because of
// the orientation of faces and lines. On the mesh of
// GridGenerator::hyper_ball_balanced(), the interpolation error of a smooth
// field has to be of the same size as the error of the L2 projection. On a
// cube mesh whose cells are rotated such that faces and lines with
// non-standard orientation appear, a field in the finite element space has to
// be reproduced exactly. Both require the values of the degrees of
// freedom to follow the sign changes of the shape functions.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
class SmoothField : public Function<dim>
{
public:
  SmoothField()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    if (component == 0)
      return std::sin(2. * p[0] + p[1]) + p[dim - 1] * p[dim - 1];
    else if (component == 1)
      return std::cos(p[0] - 3. * p[1]) * (1. + p[dim - 1]);
    else
      return std::exp(0.5 * p[0]) * std::sin(p[1] + 2. * p[2]);
  }
};



// A vector field with components in Q_k, i.e., of degree k in each
// coordinate. On affine cells, it lies in the spaces of FE_RaviartThomas,
// FE_RaviartThomasNodal, and FE_Nedelec of degree k. Each component is a
// product of generic one-dimensional polynomials, such that all face moments
// are nonzero.
template <int dim>
class TensorProductField : public Function<dim>
{
public:
  TensorProductField(const unsigned int degree)
    : Function<dim>(dim)
    , degree(degree)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    double result = 1.;
    for (unsigned int d = 0; d < dim; ++d)
      result *=
        std::pow(0.3 + 0.1 * (component + d) + (0.7 - 0.2 * d) * p[d], degree);
    return result;
  }

private:
  const unsigned int degree;
};



void
create_rotated_cube(Triangulation<3> &tria)
{
  Triangulation<3> cube;
  GridGenerator::subdivided_hyper_cube(cube, 3);
  auto [vertices, cells, subcell_data] =
    GridTools::get_coarse_mesh_description(cube);
  const auto   rots    = TestGrids::hypercube_rotations();
  unsigned int counter = 0;
  for (auto &cell : cells)
    TestGrids::rotate_cell(cell, rots[(7 * counter++ + 3) % rots.size()]);
  tria.create_triangulation(vertices, cells, SubCellData());
}



template <int dim>
double
compute_error(const Triangulation<dim> &tria,
              const Mapping<dim>       &mapping,
              const DoFHandler<dim>    &dof_handler,
              const Vector<double>     &solution,
              const Function<dim>      &function)
{
  Vector<double> cellwise_error(tria.n_active_cells());
  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    function,
                                    cellwise_error,
                                    QGauss<dim>(dof_handler.get_fe().degree +
                                                3),
                                    VectorTools::L2_norm);
  return VectorTools::compute_global_error(tria,
                                           cellwise_error,
                                           VectorTools::L2_norm);
}



// Compare the interpolation error of a smooth field to the error of the L2
// projection.
template <int dim>
void
test_smooth(const Triangulation<dim> &tria,
            const Mapping<dim>       &mapping,
            const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  const SmoothField<dim> function;
  Vector<double>         interpolated(dof_handler.n_dofs());
  VectorTools::interpolate(mapping, dof_handler, function, interpolated);

  // the hp overload of VectorTools::project() always assembles a matrix
  Vector<double> projected(dof_handler.n_dofs());
  VectorTools::project(hp::MappingCollection<dim>(mapping),
                       dof_handler,
                       constraints,
                       hp::QCollection<dim>(QGauss<dim>(fe.degree + 2)),
                       function,
                       projected);

  deallog << "  " << std::left << std::setw(42) << fe.get_name()
          << " L2 error interpolation " << std::scientific
          << std::setprecision(3)
          << compute_error(tria, mapping, dof_handler, interpolated, function)
          << ", projection "
          << compute_error(tria, mapping, dof_handler, projected, function)
          << std::defaultfloat << std::endl;
}



// Interpolate a field in the finite element space, which has to be
// reproduced exactly.
template <int dim>
void
test_exact(const Triangulation<dim> &tria,
           const Mapping<dim>       &mapping,
           const FiniteElement<dim> &fe,
           const unsigned int        degree)
{
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  const TensorProductField<dim> function(degree);
  Vector<double>                interpolated(dof_handler.n_dofs());
  VectorTools::interpolate(mapping, dof_handler, function, interpolated);

  deallog << "  " << std::left << std::setw(42) << fe.get_name()
          << " L2 error interpolation " << std::scientific
          << std::setprecision(3)
          << compute_error(tria, mapping, dof_handler, interpolated, function)
          << std::defaultfloat << std::endl;
}



int
main()
{
  initlog();

  {
    Triangulation<2> tria;
    GridGenerator::hyper_ball_balanced(tria);
    tria.refine_global(1);
    deallog << "hyper_ball_balanced, 2d, smooth field" << std::endl;
    for (unsigned int degree = 0; degree < 3; ++degree)
      {
        test_smooth(tria, MappingQ<2>(3), FE_RaviartThomas<2>(degree));
        test_smooth(tria, MappingQ<2>(3), FE_RaviartThomasNodal<2>(degree));
        test_smooth(tria,
                    MappingQ<2>(3),
                    FESystem<2>(FE_RaviartThomas<2>(degree), 1));
      }
  }
  {
    Triangulation<3> tria;
    GridGenerator::hyper_ball_balanced(tria);
    deallog << "hyper_ball_balanced, 3d, smooth field" << std::endl;
    for (unsigned int degree = 0; degree < 2; ++degree)
      {
        test_smooth(tria, MappingQ<3>(2), FE_RaviartThomas<3>(degree));
        test_smooth(tria, MappingQ<3>(2), FE_RaviartThomasNodal<3>(degree));
        test_smooth(tria,
                    MappingQ<3>(2),
                    FESystem<3>(FE_RaviartThomas<3>(degree), 1));
      }
  }
  {
    // FE_Nedelec changes signs only of dofs on lines with non-standard
    // orientation, which this mesh has
    Triangulation<3> tria;
    create_rotated_cube(tria);
    deallog << "cube with rotated cells, 3d, field in the finite element space"
            << std::endl;
    for (unsigned int degree = 0; degree < 2; ++degree)
      {
        test_exact(tria, MappingQ<3>(1), FE_RaviartThomas<3>(degree), degree);
        test_exact(tria,
                   MappingQ<3>(1),
                   FE_RaviartThomasNodal<3>(degree),
                   degree);
        test_exact(tria,
                   MappingQ<3>(1),
                   FESystem<3>(FE_RaviartThomas<3>(degree), 1),
                   degree);
        test_exact(tria, MappingQ<3>(1), FE_Nedelec<3>(degree), degree);
      }
  }
}
