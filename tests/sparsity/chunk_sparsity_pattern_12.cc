// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check the copy constructor of ChunkSparsityPattern

#include "sparsity_pattern_common.h"

int
main()
{
  initlog();

  ChunkSparsityPattern chunk1;
  ChunkSparsityPattern chunk2 = chunk1;

  deallog << "OK" << std::endl;
}
