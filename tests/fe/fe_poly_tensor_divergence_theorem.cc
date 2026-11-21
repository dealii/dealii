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
#include <array>

#include "../tests.h"


// all rotations of the unit cube as maps of the lexicographic vertex
// numbering, i.e., signed permutations of the coordinates with determinant 1
std::vector<std::array<unsigned int, 8>>
rotations()
{
  std::vector<std::array<unsigned int, 8>> result;
  std::array<unsigned int, 3>              perm = {{0, 1, 2}};
  do
    {
      int parity = 1;
      for (unsigned int a = 0; a < 3; ++a)
        for (unsigned int b = a + 1; b < 3; ++b)
          if (perm[a] > perm[b])
            parity = -parity;
      for (unsigned int signs = 0; signs < 8; ++signs)
        {
          int det = parity;
          for (unsigned int d = 0; d < 3; ++d)
            if (signs & (1u << d))
              det = -det;
          if (det < 0)
            continue;
          std::array<unsigned int, 8> map;
          for (unsigned int v = 0; v < 8; ++v)
            {
              unsigned int w = 0;
              for (unsigned int d = 0; d < 3; ++d)
                {
                  unsigned int bit = (v >> perm[d]) & 1;
                  if (signs & (1u << d))
                    bit = 1 - bit;
                  w |= bit << d;
                }
              map[v] = w;
            }
          result.push_back(map);
        }
    }
  while (std::next_permutation(perm.begin(), perm.end()));
  return result;
}



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
  const auto   rots    = rotations();
  unsigned int counter = 0;
  for (auto &cell : cells)
    {
      const auto &rotation = rots[(7 * counter++ + 3) % rots.size()];
      std::array<unsigned int, 8> original;
      for (unsigned int v = 0; v < 8; ++v)
        original[v] = cell.vertices[v];
      for (unsigned int v = 0; v < 8; ++v)
        cell.vertices[rotation[v]] = original[v];
    }
  tria.create_triangulation(vertices, cells, SubCellData());
}



template <int dim>
void
test(const Triangulation<dim> &tria,
     const FiniteElement<dim> &fe,
     const bool                check_hessians)
{
  // a linear mapping, such that the cells are not parallelograms but the
  // quadrature formulas below integrate accurately
  const MappingQ<dim> mapping(1);
  DoFHandler<dim>     dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  FEValues<dim>     fe_values(mapping,
                              fe,
                              QGauss<dim>(fe.degree + 8),
                              update_gradients | update_JxW_values |
                                (check_hessians ? update_hessians :
                                                  update_default));
  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   QGauss<dim - 1>(fe.degree + 8),
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
