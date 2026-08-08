// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test the support points of FE_WedgeP for consistency by checking if they
// correctly on a vertex, edge, face or in the volume

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"


void
test(const unsigned int degree)
{
  constexpr int dim = 3;
  const double  tol = 1e-12;

  deallog << "Support points of degree " << degree << std::endl;
  const auto fe             = FE_WedgeP<dim>(degree);
  const auto support_points = fe.get_unit_support_points();
  const auto reference_cell = fe.reference_cell();

  const unsigned int n_dof_per_line = degree - 1;
  const unsigned int n_dof_per_tri  = (degree - 2) * (degree - 1) / 2;
  const unsigned int n_dof_per_quad = n_dof_per_line * n_dof_per_line;
  const unsigned int n_dof_total =
    (degree + 1) * (degree + 1) * (degree + 2) / 2;

  unsigned int counter = 0;
  for (const auto &p : support_points)
    {
      const double x = p[0];
      const double y = p[1];
      const double z = p[2];

      unsigned int n_zeros = 0;
      for (unsigned int d = 0; d < dim; ++d)
        if (std::abs(p[d]) < tol)
          ++n_zeros;

      if (counter < 3)
        {
          // check vertices
          // on the bottom two coordinates have to be 0
          Assert(n_zeros > 1, ExcInternalError());
          Assert(std::abs(z) < tol, ExcInternalError());
        }
      else if (counter < 6)
        {
          // check vertices
          // one coordinate has to be 0
          Assert(n_zeros > 0, ExcInternalError());
          Assert(std::abs(z - 1.0) < tol, ExcInternalError());
        }
      else if (counter < 6 + 1 * n_dof_per_line)
        {
          // first line
          Assert(n_zeros > 1, ExcInternalError());
          Assert(std::abs(y) < tol, ExcInternalError());
          Assert(std::abs(z) < tol, ExcInternalError());
        }
      else if (counter < 6 + 2 * n_dof_per_line)
        {
          // second line
          Assert(n_zeros > 0, ExcInternalError());
          Assert(std::abs(x + y - 1.0) < tol, ExcInternalError());
          Assert(std::abs(z) < tol, ExcInternalError());
        }
      else if (counter < 6 + 3 * n_dof_per_line)
        {
          // third line
          Assert(n_zeros > 0, ExcInternalError());
          Assert(std::abs(x) < tol, ExcInternalError());
          Assert(std::abs(z) < tol, ExcInternalError());
        }
      else if (counter < 6 + 4 * n_dof_per_line)
        {
          // 4th line
          Assert(n_zeros > 0, ExcInternalError());
          Assert(std::abs(y) < tol, ExcInternalError());
          Assert(std::abs(z - 1.0) < tol, ExcInternalError());
        }
      else if (counter < 6 + 5 * n_dof_per_line)
        {
          // 5th line
          Assert(std::abs(x + y - 1.0) < tol, ExcInternalError());
          Assert(std::abs(z - 1.0) < tol, ExcInternalError());
        }
      else if (counter < 6 + 6 * n_dof_per_line)
        {
          // 6th line
          Assert(n_zeros > 0, ExcInternalError());
          Assert(std::abs(x) < tol, ExcInternalError());
          Assert(std::abs(z - 1.0) < tol, ExcInternalError());
        }
      else if (counter < 6 + 7 * n_dof_per_line)
        {
          // 7th line
          Assert(std::abs(x) < tol, ExcInternalError());
          Assert(std::abs(y) < tol, ExcInternalError());
        }
      else if (counter < 6 + 8 * n_dof_per_line)
        {
          // 8th line
          Assert(std::abs(x - 1.0) < tol, ExcInternalError());
          Assert(std::abs(y) < tol, ExcInternalError());
        }
      else if (counter < 6 + 9 * n_dof_per_line)
        {
          // 9th line
          Assert(std::abs(x) < tol, ExcInternalError());
          Assert(std::abs(y - 1.0) < tol, ExcInternalError());
        }
      else if (counter < 6 + 9 * n_dof_per_line + 1 * n_dof_per_tri)
        {
          // first face
          Assert(std::abs(z) < tol, ExcInternalError());
          Assert(n_zeros == 1, ExcInternalError());
        }
      else if (counter < 6 + 9 * n_dof_per_line + 2 * n_dof_per_tri)
        {
          // second face
          Assert(std::abs(z - 1.0) < tol, ExcInternalError());
          Assert(n_zeros == 0, ExcInternalError());
        }
      else if (counter <
               6 + 9 * n_dof_per_line + 2 * n_dof_per_tri + 1 * n_dof_per_quad)
        {
          // third face
          Assert(std::abs(y) < tol, ExcInternalError());
          Assert(n_zeros == 1, ExcInternalError());
        }
      else if (counter <
               6 + 9 * n_dof_per_line + 2 * n_dof_per_tri + 2 * n_dof_per_quad)
        {
          // 4th face
          Assert(std::abs(x + y - 1.0) < tol, ExcInternalError());
          Assert(n_zeros == 0, ExcInternalError());
        }
      else if (counter <
               6 + 9 * n_dof_per_line + 2 * n_dof_per_tri + 3 * n_dof_per_quad)
        {
          // 5th face
          Assert(std::abs(x) < tol, ExcInternalError());
          Assert(n_zeros == 1, ExcInternalError());
        }
      else if (counter < n_dof_total)
        {
          // volume
          Assert(n_zeros == 0, ExcInternalError());
        }
      else
        DEAL_II_ASSERT_UNREACHABLE();

      ++counter;
    }

  deallog << "All nodes correct" << std::endl;

  deallog << std::endl;
}


int
main()
{
  initlog();

  deallog.push("3D");
  for (unsigned int i = 1; i < 3; ++i)
    test(i);
  deallog.pop();
}
