
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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>


template <int dim, typename NumberType>
Tensor<2, dim, NumberType>
make_tensor(const NumberType &val)
{
  Tensor<2, dim, NumberType> out;
  for (unsigned int i = 0; i < dim; ++i)
    out[i][i] = 1.0;

  for (unsigned int i = 0; i < out.n_independent_components; ++i)
    out[out.unrolled_to_component_indices(i)] += i + val;
  return out;
}


template <int dim, typename NumberType>
SymmetricTensor<2, dim, NumberType>
make_symm_tensor(const NumberType &val)
{
  SymmetricTensor<2, dim, NumberType> out;
  for (unsigned int i = 0; i < dim; ++i)
    out[i][i] = 1.0;

  for (unsigned int i = 0; i < out.n_independent_components; ++i)
    out[out.unrolled_to_component_indices(i)] += i + val;
  return out;
}


template <typename Stream, typename NumberType>
void
print(Stream &stream, const std::string &name, const NumberType &val)
{
  stream << name << ": " << val << std::endl;
}


template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<0, dim, NumberType> &val)
{
  stream << name << ": " << val << std::endl;
}


template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<1, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    stream << name << "[" << i << "]: " << t[i] << std::endl;
}


template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<2, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      stream << name << "[" << i << "][" << j << "]: " << t[i][j] << std::endl;
}


template <typename Stream, int dim, typename NumberType>
void
print(Stream &                                   stream,
      const std::string &                        name,
      const SymmetricTensor<2, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      stream << name << "[" << i << "][" << j << "]: " << t[i][j] << std::endl;
}


template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<3, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        stream << name << "[" << i << "][" << j << "][" << k
               << "]: " << t[i][j][k] << std::endl;
}


template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<4, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          stream << name << "[" << i << "][" << j << "][" << k << "][" << l
                 << "]: " << t[i][j][k][l] << std::endl;
}


template <typename Stream, int dim, typename NumberType>
void
print(Stream &                                   stream,
      const std::string &                        name,
      const SymmetricTensor<4, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          stream << name << "[" << i << "][" << j << "][" << k << "][" << l
                 << "]: " << t[i][j][k][l] << std::endl;
}
