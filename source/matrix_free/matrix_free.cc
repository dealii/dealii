// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2016 by the deal.II authors
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


#include <deal.II/matrix_free/matrix_free.templates.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <iostream>

DEAL_II_NAMESPACE_OPEN

#include "matrix_free.inst"

template struct internal::MatrixFreeFunctions::ShapeInfo<double>;
template struct internal::MatrixFreeFunctions::ShapeInfo<float>;

DEAL_II_NAMESPACE_CLOSE
