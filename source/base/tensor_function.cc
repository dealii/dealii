// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2015 by the deal.II authors
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

#include <deal.II/base/tensor_function.templates.h>

DEAL_II_NAMESPACE_OPEN

// Explicit instantiations:

template class TensorFunction<1,1>;
template class TensorFunction<2,1>;
template class TensorFunction<3,1>;
template class TensorFunction<4,1>;
template class TensorFunction<1,2>;
template class TensorFunction<2,2>;
template class TensorFunction<3,2>;
template class TensorFunction<4,2>;
template class TensorFunction<1,3>;
template class TensorFunction<2,3>;
template class TensorFunction<3,3>;
template class TensorFunction<4,3>;

template class ConstantTensorFunction<1,1>;
template class ConstantTensorFunction<2,1>;
template class ConstantTensorFunction<3,1>;
template class ConstantTensorFunction<4,1>;
template class ConstantTensorFunction<1,2>;
template class ConstantTensorFunction<2,2>;
template class ConstantTensorFunction<3,2>;
template class ConstantTensorFunction<4,2>;
template class ConstantTensorFunction<1,3>;
template class ConstantTensorFunction<2,3>;
template class ConstantTensorFunction<3,3>;
template class ConstantTensorFunction<4,3>;

template class ZeroTensorFunction<1,1>;
template class ZeroTensorFunction<2,1>;
template class ZeroTensorFunction<3,1>;
template class ZeroTensorFunction<4,1>;
template class ZeroTensorFunction<1,2>;
template class ZeroTensorFunction<2,2>;
template class ZeroTensorFunction<3,2>;
template class ZeroTensorFunction<4,2>;
template class ZeroTensorFunction<1,3>;
template class ZeroTensorFunction<2,3>;
template class ZeroTensorFunction<3,3>;
template class ZeroTensorFunction<4,3>;


DEAL_II_NAMESPACE_CLOSE
