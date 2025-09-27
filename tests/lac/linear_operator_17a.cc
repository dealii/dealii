// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test the LinearOperator template on a trivial vector implementation
// :: RightVector -> LeftVector

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/linear_operator.h>

#include "../tests.h"



int
main()
{
  initlog();

  const std::function<void(Vector<double> &, bool)> reinit_vector =
    [](Vector<double> &v, bool omit_zeroing_entries) {
      v.reinit(3, omit_zeroing_entries);
    };

  auto       id     = identity_operator(reinit_vector);
  const auto filter = mean_value_filter(id);

  const auto block_filter = block_diagonal_operator<3>(filter);

  BlockVector<double> block_vector(3, 3);
  block_vector.block(0)[0] = 1;
  block_vector.block(0)[1] = 2;
  block_vector.block(0)[2] = 3;
  block_vector.block(1)[0] = 4;
  block_vector.block(1)[1] = 6;
  block_vector.block(1)[2] = 8;
  block_vector.block(2)[0] = 9;
  block_vector.block(2)[1] = 12;
  block_vector.block(2)[2] = 15;

  deallog << block_vector.block(0) << std::endl;
  deallog << block_vector.block(1) << std::endl;
  deallog << block_vector.block(2) << std::endl;

  BlockVector<double> block_vector_2(3, 3);
  block_filter.vmult(block_vector_2, block_vector);

  deallog << block_vector_2.block(0) << std::endl;
  deallog << block_vector_2.block(1) << std::endl;
  deallog << block_vector_2.block(2) << std::endl;
}
