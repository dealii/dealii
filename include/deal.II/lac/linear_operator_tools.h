// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_linear_operator_tools_h
#define dealii_linear_operator_tools_h

// Many usage cases lead to a combination of LinearOperator and
// PackagedOperation. To ease the pain of reading compilation errors, just
// include all headers we ever need to use LO and friends in one place:

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/constrained_linear_operator.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/schur_complement.h>
#include <deal.II/lac/trilinos_linear_operator.h>

#endif
