// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Compress/decompress ConstraintKinds.

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/evaluation_kernels_hanging_nodes.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

int
main()
{
  initlog();

  using namespace dealii::internal::MatrixFreeFunctions;

  const unsigned int dim = 3;

  const auto check = [&](const auto subcell, const auto face, const auto edge) {
    const auto kind =
      static_cast<ConstraintKinds>(subcell + (face << 3) + (edge << 6));

    Assert(kind == decompress(compress(kind, dim), dim), ExcInternalError());
  };

  check(0, 0, 0);

  for (unsigned int subcell = 0; subcell < 8; ++subcell)
    {
      for (unsigned int face = 1; face < 8; ++face)
        check(subcell, face, 0);

      for (unsigned int edge = 1; edge < 8; ++edge)
        check(subcell, 0, edge);

      for (unsigned int d = 0; d < dim; ++d)
        check(subcell, 1 << d, 1 << d);
    }

  deallog << "OK!" << std::endl;
}
