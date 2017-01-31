// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


#ifndef dealii_tensor_classes_fwd_h
#define dealii_tensor_classes_fwd_h

#include <deal.II/base/config.h>

// C++ does not permit repeated forward declarations of template classes with
// default template arguments. Hence, collect the four most common
// Tensor-like classes here so that we can just include this file (with the
// <code>#ifndef</code> guards) to avoid this problem.
DEAL_II_NAMESPACE_OPEN

template <int spacedim, typename Number = double>
class Point;
template <int rank, int dim, typename Number = double>
class Tensor;
template <int rank, int dim, typename Number = double>
class SymmetricTensor;
template <int order, int dim, int spacedim, typename Number = double>
class DerivativeForm;

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_tensor_classes_fwd_h
