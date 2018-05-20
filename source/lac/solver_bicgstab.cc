// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#include <deal.II/lac/solver_bicgstab.h>

DEAL_II_NAMESPACE_OPEN

internal::SolverBicgstabData::SolverBicgstabData()
  : alpha(0.),
    beta(0.),
    omega(0.),
    rho(0.),
    rhobar(0.),
    step(numbers::invalid_unsigned_int),
    res(numbers::signaling_nan<double>())
{}

DEAL_II_NAMESPACE_CLOSE
