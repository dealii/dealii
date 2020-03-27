// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


// Test that rolling to, and unrolling from, Kelvin notation.

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/physics/notation.h>

#include "../tests.h"


using namespace dealii::Physics;

template <int dim, typename Number>
void initialize(Tensor<1, dim, Number> &x)
{
  for (unsigned int i = 0; i < x.n_independent_components; ++i)
    x[i] = i + 1;
}

template <int dim, typename Number>
void initialize(Tensor<2, dim, Number> &x)
{
  unsigned int c = 1;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        x[i][j] = c++;
      }
}

template <int dim, typename Number>
void initialize(SymmetricTensor<2, dim, Number> &x)
{
  unsigned int c = 1;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      {
        x[i][j] = c++;
      }
}

template <int dim, typename Number>
void initialize(Tensor<3, dim, Number> &x)
{
  unsigned int c = 1;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        {
          x[i][j][k] = c++;
        }
}

// A specialised constructor mimicking the construction
// of a rank-3 tensor with two symmetric components
template <int dim, typename Number>
void initialize(Tensor<3, dim, Number> &x,
                const bool              left_components_are_symmetric)
{
  Tensor<1, dim, Number> v;
  initialize(v);
  SymmetricTensor<2, dim, Number> st;
  initialize(st);
  const Tensor<2, dim, Number> t(st);
  if (left_components_are_symmetric == true)
    x = outer_product(t, v);
  else
    x = outer_product(v, t);
}

template <int dim, typename Number>
void initialize(Tensor<4, dim, Number> &x)
{
  unsigned int c = 1;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          {
            x[i][j][k][l] = c++;
          }
}

template <int dim, typename Number>
void initialize(SymmetricTensor<4, dim, Number> &x)
{
  unsigned int c = 1;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          {
            x[i][j][k][l] = c++;
          }
}

template <int dim>
void
test_scalars()
{
  const double A = 5;

  const Vector<double>     vA = Notation::Kelvin::to_vector(A);
  const FullMatrix<double> mA = Notation::Kelvin::to_matrix(A);

  typedef typename std::decay<decltype(A)>::type InpType;
  const auto vA_conv = Notation::Kelvin::to_tensor<InpType>(vA);
  const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

  Assert(std::abs(vA_conv - A) < 1e-12,
         ExcMessage("Different result for vector conversion"));
  Assert(std::abs(mA_conv - A) < 1e-12,
         ExcMessage("Different result for matrix conversion"));
}

template <int dim>
void
test_rank_0_tensors()
{
  const Tensor<0, dim, double> A = 5;

  const Vector<double>     vA = Notation::Kelvin::to_vector(A);
  const FullMatrix<double> mA = Notation::Kelvin::to_matrix(A);

  typedef typename std::decay<decltype(A)>::type InpType;
  const auto vA_conv = Notation::Kelvin::to_tensor<InpType>(vA);
  const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

  Assert(std::abs(vA_conv - A) < 1e-12,
         ExcMessage("Different result for vector conversion"));
  Assert(std::abs(mA_conv - A) < 1e-12,
         ExcMessage("Different result for matrix conversion"));
}

template <int dim>
void
test_rank_1_tensors()
{
  Tensor<1, dim, double> A;
  initialize(A);

  const Vector<double> vA = Notation::Kelvin::to_vector(A);

  typedef typename std::decay<decltype(A)>::type InpType;
  const auto vA_conv = Notation::Kelvin::to_tensor<InpType>(vA);

  Assert((vA_conv - A).norm() < 1e-12,
         ExcMessage("Different result for vector conversion"));
}

template <int dim>
void
test_rank_2_tensors()
{
  // Non-symmetric tensor
  {
    Tensor<2, dim, double> A;
    initialize(A);

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mA = Notation::Kelvin::to_matrix(A);

    typedef typename std::decay<decltype(A)>::type InpType;
    const auto vA_conv = Notation::Kelvin::to_tensor<InpType>(vA);
    const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

    Assert((vA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for vector conversion"));
    Assert((mA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for matrix conversion"));
  }

  // Symmetric tensor
  {
    SymmetricTensor<2, dim, double> A;
    initialize(A);

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mA = Notation::Kelvin::to_matrix(A);

    typedef typename std::decay<decltype(A)>::type InpType;
    const auto vA_conv = Notation::Kelvin::to_tensor<InpType>(vA);
    const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

    Assert((vA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for vector conversion"));
    Assert((mA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for matrix conversion"));
  }
}

template <int dim>
void
test_rank_3_tensors()
{
  // Non-symmetric tensor: Version 1
  {
    Tensor<3, dim, double> A;
    initialize(A);

    const FullMatrix<double> mA =
      Notation::Kelvin::to_matrix<dim, Tensor<1, dim>, Tensor<2, dim>>(A);

    typedef typename std::decay<decltype(A)>::type InpType;
    const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

    Assert((mA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for matrix conversion"));
  }

  // Non-symmetric tensor: Version 2
  {
    Tensor<3, dim, double> A;
    initialize(A);

    const FullMatrix<double> mA =
      Notation::Kelvin::to_matrix<dim, Tensor<2, dim>, Tensor<1, dim>>(A);

    typedef typename std::decay<decltype(A)>::type InpType;
    const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

    Assert((mA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for matrix conversion"));
  }

  // Symmetric tensor: Version 1
  {
    Tensor<3, dim, double> A;
    initialize(A, true); // Specialised constructor

    const FullMatrix<double> mA =
      Notation::Kelvin::to_matrix<dim, SymmetricTensor<2, dim>, Tensor<1, dim>>(
        A);

    typedef typename std::decay<decltype(A)>::type InpType;
    const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

    Assert((mA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for matrix conversion"));
  }

  // Symmetric tensor: Version 2
  {
    Tensor<3, dim, double> A;
    initialize(A, false); // Specialised constructor

    const FullMatrix<double> mA =
      Notation::Kelvin::to_matrix<dim, Tensor<1, dim>, SymmetricTensor<2, dim>>(
        A);

    typedef typename std::decay<decltype(A)>::type InpType;
    const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

    Assert((mA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for matrix conversion"));
  }
}

template <int dim>
void
test_rank_4_tensors()
{
  // Non-symmetric tensor
  {
    Tensor<4, dim, double> A;
    initialize(A);

    const FullMatrix<double> mA = Notation::Kelvin::to_matrix(A);

    typedef typename std::decay<decltype(A)>::type InpType;
    const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

    Assert((mA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for matrix conversion"));
  }

  // Symmetric tensor
  {
    SymmetricTensor<4, dim, double> A;
    initialize(A);

    const FullMatrix<double> mA = Notation::Kelvin::to_matrix(A);

    typedef typename std::decay<decltype(A)>::type InpType;
    const auto mA_conv = Notation::Kelvin::to_tensor<InpType>(mA);

    Assert((mA_conv - A).norm() < 1e-12,
           ExcMessage("Different result for matrix conversion"));
  }
}

template <int dim>
void
test_tensors()
{
  test_scalars<dim>();
  test_rank_0_tensors<dim>();
  test_rank_1_tensors<dim>();
  test_rank_2_tensors<dim>();
  test_rank_3_tensors<dim>();
  test_rank_4_tensors<dim>();
}

int
main()
{
  initlog();

  test_tensors<2>();
  test_tensors<3>();

  deallog << "OK" << std::endl;
}
