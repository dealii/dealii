// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2023 by the deal.II authors
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


#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_internal.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.templates.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

// explicit instantiation
#include "mg_transfer_matrix_free.inst"


DEAL_II_NAMESPACE_CLOSE
