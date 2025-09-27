// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
//
// Test a tolerance parameter to find the periodic nodes of a triangulation
// using the collect_periodic_faces function.
// Useful in case of CAD geometry that is not precise enough. The default
// value is at 1e-10 (absolute value). Sequential version.

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <vector>

#include "../tests.h"

template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;
  Triangulation<dim> tria;

  std::vector<typename GridTools::PeriodicFacePair<
    typename Triangulation<dim>::cell_iterator>>
    periodicity_vector;

  unsigned int              num_refinements = 1 << 4;
  Point<dim>                p1;
  Point<dim>                p2;
  std::vector<unsigned int> repetitions(dim);

  for (unsigned int i = 0; i < dim; ++i)
    {
      p1[i]          = 0.;
      p2[i]          = 1.;
      repetitions[i] = num_refinements;
    }
  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2, true);

  // Here we are only interested in using the tolerance parameter of the
  // collect_periodic_faces function, the other default parameter stay as
  // default one.
  const Tensor<1, dim>     defaultOffset = dealii::Tensor<1, dim>();
  const FullMatrix<double> matrix        = FullMatrix<double>();
  const double             tolerance     = 1e-7;
  const unsigned int       direction     = 0; // x-direction

  // We distort the mesh to test the tolerance parameter
  GridTools::distort_random(tolerance / 10, tria, /* keep boundary */ false);

  // Collect periodic faces in the x-direction
  // activate debug to see it fail
  GridTools::collect_periodic_faces(tria,
                                    0, // boundary Id 1
                                    1, // boundary Id 2
                                    direction,
                                    periodicity_vector,
                                    defaultOffset,
                                    matrix,
                                    tolerance);
  // Check the size of the matched_pairs vector
  deallog << periodicity_vector.size() << std::endl;
}

int
main(int argc, char *argv[])
{
  initlog();

  test<2>();
  test<3>();

  return 0;
}
