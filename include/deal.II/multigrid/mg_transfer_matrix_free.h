// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mg_transfer_matrix_free_h
#define dealii_mg_transfer_matrix_free_h

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>


DEAL_II_NAMESPACE_OPEN


template <int dim, typename Number>
using MGTransferMatrixFree = MGTransferMF<dim, Number>;

template <int dim, typename Number>
using MGTransferBlockMatrixFree = MGTransferBlockMF<dim, Number>;


DEAL_II_NAMESPACE_CLOSE

#endif
