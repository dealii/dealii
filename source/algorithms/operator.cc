// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2014 by the deal.II authors
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


#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/algorithms/operator.templates.h>
#include <deal.II/algorithms/newton.templates.h>
#include <deal.II/algorithms/theta_timestepping.templates.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>


DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  OperatorBase::~OperatorBase()
  {}

  void OperatorBase::notify(const Event &e)
  {
    notifications += e;
  }

  void
  OperatorBase::clear_events ()
  {
    notifications.clear();
  }


#include "operator.inst"
}

DEAL_II_NAMESPACE_CLOSE
