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

// Like fe_rt_gradient_divergence_theorem and fe_rt_hessian_divergence_theorem,
// check the shape gradients and Hessians of FE_PolyTensor elements by the
// divergence theorem: the integral of the gradient (Hessian) of each shape
// function over a cell has to match the integral of the shape function
// (its gradient) times the normal vector over the boundary of the cell. The
// cells are not parallelograms and rotated against each other, such that
// FEValues changes the sign of some shape functions: the mesh of
// GridGenerator::hyper_ball_balanced() in 2d and a deformed cube mesh in 3d
// whose cells are rotated such that all face and line orientations appear.
// The terms of the gradients and Hessians that come from the derivatives of
// the Jacobian did not get the sign change. The Hessians are only checked for
// FE_Nedelec: those of the H(div) elements are not correct on cells that are
// not parallelograms, also without a sign change.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <algorithm>

#include "../tests.h"

#include "../test_grids.h"


void
create_rotated_cube(Triangulation<3> &tria)
{
  Triangulation<3> cube;
  GridGenerator::subdivided_hyper_cube(cube, 3);
  GridTools::transform(
    [](const Point<3> &p) {
      Point<3> q = p;
      q[0] += 0.08 * std::sin(3. * p[1] + 1.);
      q[1] += 0.06 * std::cos(2. * p[0] + p[2]);
      q[2] += 0.05 * std::sin(2. * p[0] + 3. * p[1]);
      return q;
    },
    cube);
  auto [vertices, cells, subcell_data] =
    GridTools::get_coarse_mesh_description(cube);
  const auto   rots    = TestGrids::hypercube_rotations();
  unsigned int counter = 0;
  for (auto &cell : cells)
    TestGrids::rotate_cell(cell, rots[(7 * counter++ + 3) % rots.size()]);
  tria.create_triangulation(vertices, cells, SubCellData());
}



template <int dim>
void
test(const Triangulation<dim> &tria,
     const FiniteElement<dim> &fe,
     const bool                check_hessians)
{
  // a linear mapping, such that the cells are not parallelograms. The
  // integrands are then not polynomials: the quadrature formulas below use
  // one point more than needed to meet the tolerance below, the 2d mesh needs
  // the most points for FE_RaviartThomasNodal of degree 2
  const MappingQ<dim> mapping(1);
  const unsigned int  n_points = fe.degree + (dim == 2 ? 8 : 3);
  DoFHandler<dim>     dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  FEValues<dim>     fe_values(mapping,
                              fe,
                              QGauss<dim>(n_points),
                              update_gradients | update_JxW_values |
                                (check_hessians ? update_hessians :
                                                  update_default));
  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   QGauss<dim - 1>(n_points),
                                   update_values | update_normal_vectors |
                                     update_JxW_values |
                                     (check_hessians ? update_gradients :
                                                       update_default));

  const FEValuesExtractors::Vector vector(0);
  std::vector<Tensor<2, dim>>      boundary_gradients(fe.n_dofs_per_cell());
  std::vector<Tensor<3, dim>>      boundary_hessians(fe.n_dofs_per_cell());
  unsigned int n_failed_cells_gradients = 0, n_failed_cells_hessians = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      std::fill(boundary_gradients.begin(),
                boundary_gradients.end(),
                Tensor<2, dim>());
      std::fill(boundary_hessians.begin(),
                boundary_hessians.end(),
                Tensor<3, dim>());
      for (const unsigned int face : cell->face_indices())
        {
          fe_face_values.reinit(cell, face);
          for (const unsigned int i : fe_face_values.dof_indices())
            for (const unsigned int q :
                 fe_face_values.quadrature_point_indices())
              {
                boundary_gradients[i] +=
                  outer_product(fe_face_values[vector].value(i, q),
                                fe_face_values.normal_vector(q)) *
                  fe_face_values.JxW(q);
                if (check_hessians)
                  boundary_hessians[i] +=
                    outer_product(fe_face_values[vector].gradient(i, q),
                                  fe_face_values.normal_vector(q)) *
                    fe_face_values.JxW(q);
              }
        }

      fe_values.reinit(cell);
      bool gradients_ok = true, hessians_ok = true;
      for (const unsigned int i : fe_values.dof_indices())
        {
          Tensor<2, dim> bulk_gradient;
          Tensor<3, dim> bulk_hessian;
          for (const unsigned int q : fe_values.quadrature_point_indices())
            {
              bulk_gradient +=
                fe_values[vector].gradient(i, q) * fe_values.JxW(q);
              if (check_hessians)
                bulk_hessian +=
                  fe_values[vector].hessian(i, q) * fe_values.JxW(q);
            }
          if ((bulk_gradient - boundary_gradients[i]).norm() >
              1e-8 * std::max(boundary_gradients[i].norm(), 1e-6))
            gradients_ok = false;
          if ((bulk_hessian - boundary_hessians[i]).norm() >
              1e-8 * std::max(boundary_hessians[i].norm(), 1e-6))
            hessians_ok = false;
        }
      if (!gradients_ok)
        ++n_failed_cells_gradients;
      if (!hessians_ok)
        ++n_failed_cells_hessians;
    }

  deallog << "  " << fe.get_name() << ": integrals of the gradients do not "
          << "match on " << n_failed_cells_gradients << " cells";
  if (check_hessians)
    deallog << ", of the Hessians on " << n_failed_cells_hessians << " cells";
  deallog << std::endl;
}



int
main()
{
  initlog();

  {
    Triangulation<2> tria;
    GridGenerator::hyper_ball_balanced(tria);
    deallog << "hyper_ball_balanced, 2d" << std::endl;
    for (unsigned int degree = 0; degree < 3; ++degree)
      {
        test<2>(tria, FE_RaviartThomas<2>(degree), false);
        test<2>(tria, FE_RaviartThomasNodal<2>(degree), false);
        test<2>(tria, FE_Nedelec<2>(degree), true);
      }
  }
  {
    Triangulation<3> tria;
    create_rotated_cube(tria);
    deallog << "deformed cube with rotated cells, 3d" << std::endl;
    for (unsigned int degree = 0; degree < 2; ++degree)
      {
        test<3>(tria, FE_RaviartThomas<3>(degree), false);
        test<3>(tria, FE_RaviartThomasNodal<3>(degree), false);
        test<3>(tria, FE_Nedelec<3>(degree), true);
      }
  }
}
