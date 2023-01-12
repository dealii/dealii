// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Test consistency of ReferenceCell::volume() and
// ReferenceCell::contains_point() by computing a Monte Carlo integral
// of the volume based on the characteristic function implied by
// contains_point().

#include <deal.II/grid/reference_cell.h>

#include <random>

#include "../tests.h"


using namespace dealii;

template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  unsigned int       n_samples_inside = 0;
  const unsigned int n_samples        = 200000;

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
