// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2015 by the deal.II authors
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

#ifndef dealii__compressed_simple_sparsity_pattern_h
#define dealii__compressed_simple_sparsity_pattern_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <vector>
#include <algorithm>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

template <typename number> class SparseMatrix;


/*! @addtogroup Sparsity
 *@{
 */

/**
 * @deprecated Use DynamicSparsityPattern instead.
 */
typedef DynamicSparsityPattern CompressedSimpleSparsityPattern DEAL_II_DEPRECATED;

/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
