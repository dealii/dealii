// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#include <deal.II/base/function.templates.h>

DEAL_II_NAMESPACE_OPEN

// explicit instantiations

template class Function<1>;
template class ZeroFunction<1>;
template class ConstantFunction<1>;
template class ComponentSelectFunction<1>;
template class ScalarFunctionFromFunctionObject<1>;
template class VectorFunctionFromScalarFunctionObject<1>;
template class VectorFunctionFromTensorFunction<1>;

template class Function<2>;
template class ZeroFunction<2>;
template class ConstantFunction<2>;
template class ComponentSelectFunction<2>;
template class ScalarFunctionFromFunctionObject<2>;
template class VectorFunctionFromScalarFunctionObject<2>;
template class VectorFunctionFromTensorFunction<2>;

template class Function<3>;
template class ZeroFunction<3>;
template class ConstantFunction<3>;
template class ComponentSelectFunction<3>;
template class ScalarFunctionFromFunctionObject<3>;
template class VectorFunctionFromScalarFunctionObject<3>;
template class VectorFunctionFromTensorFunction<3>;

DEAL_II_NAMESPACE_CLOSE
