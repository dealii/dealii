// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test for constexpr EvaluationFlags

#include <deal.II/matrix_free/evaluation_flags.h>

#include "../tests.h"


int
main()
{
  initlog();

  constexpr auto flags_values_gradients =
    EvaluationFlags::values | EvaluationFlags::gradients;
  constexpr auto flags_values_hessians =
    EvaluationFlags::values | EvaluationFlags::hessians;

  static_assert(flags_values_gradients == 0x3);
  static_assert(flags_values_hessians == 0x5);
  static_assert((flags_values_gradients & flags_values_hessians) ==
                EvaluationFlags::values);
  static_assert((flags_values_gradients & EvaluationFlags::hessians) ==
                EvaluationFlags::nothing);

  deallog << "OK" << std::endl;
}
