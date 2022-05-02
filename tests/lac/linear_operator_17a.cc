// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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
