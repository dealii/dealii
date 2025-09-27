// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/algorithms/newton.templates.h>
#include <deal.II/algorithms/operator.templates.h>
#include <deal.II/algorithms/theta_timestepping.templates.h>

#include <deal.II/base/logstream.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>


DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  void
  OperatorBase::notify(const Event &e)
  {
    notifications += e;
  }



  void
  OperatorBase::clear_events()
  {
    notifications.clear();
  }

#include "algorithms/operator.inst"
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE
