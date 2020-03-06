// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#include <deal.II/base/array_view.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>

DEAL_II_NAMESPACE_OPEN

#include "mpi_noncontiguous_partitioner.inst"

DEAL_II_NAMESPACE_CLOSE
