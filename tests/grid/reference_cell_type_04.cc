// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test consistency of ReferenceCell::volume() and
// ReferenceCell::contains_point() by computing a Monte Carlo integral
// of the volume based on the characteristic function implied by
// contains_point().

#include <deal.II/grid/reference_cell.h>

#include <random>

#include "../tests.h"



template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  unsigned int       n_samples_inside = 0;
  const unsigned int n_samples        = 200000;

  // sanity check: does the reference cell contain its own nodes?
  for (const unsigned int vertex_no : reference_cell.vertex_indices())
    AssertThrow(reference_cell.contains_point(
                  reference_cell.template vertex<dim>(vertex_no)),
                ExcInternalError());

  for (unsigned int n = 0; n < n_samples; ++n)
    {
      // Choose a random point in the box [-1,1]^d that contains all
      // of our reference cells:
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = random_value<double>(-1, 1);

      if (reference_cell.contains_point(p))
        ++n_samples_inside;
    }

  // The approximate volume of the reference cell is then the fraction
  // of points that were found to be within the reference cell times
  // the volume of the box within which we have sampled:
  const double volume = 1. * n_samples_inside / n_samples * std::pow(2.0, dim);

  deallog << "ReferenceCell: " << reference_cell.to_string() << std::endl;
  deallog << "  self-reported volume = " << reference_cell.volume()
          << std::endl;
  deallog << "  computed approximate volume = " << volume << std::endl;

  Assert(std::fabs(volume - reference_cell.volume()) < 1e-2,
         ExcInternalError());
}

int
main()
{
  initlog();

  {
    deallog.push("1D");
    test<1>(ReferenceCells::Line);
    deallog.pop();
  }

  {
    deallog.push("2D");
    test<2>(ReferenceCells::Quadrilateral);
    test<2>(ReferenceCells::Triangle);
    deallog.pop();
  }

  {
    deallog.push("3D");
    test<3>(ReferenceCells::Tetrahedron);
    test<3>(ReferenceCells::Pyramid);
    test<3>(ReferenceCells::Wedge);
    test<3>(ReferenceCells::Hexahedron);
    deallog.pop();
  }
}
