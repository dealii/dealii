// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_linear_operator_tools_h
#define dealii_linear_operator_tools_h

// Many usage cases lead to a combination of LinearOperator and
// PackagedOperation. To ease the pain of reading compilation errors, just
// include all headers we ever need to use LO and friends in one place:

#include <deal.II/base/config.h>

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/constrained_linear_operator.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/schur_complement.h>
#include <deal.II/lac/trilinos_linear_operator.h>

#endif
