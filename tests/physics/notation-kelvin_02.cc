// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


// Test that Kelvin notation works as expected. The following properties
// should hold:
// - B = C:A  --> kelvin(B) = kelvin(C)*kelvin(A)
// - A:A == kelvin(A)*kelvin(A)

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
  const double C = 12;
  const double B = C * A;

  const Vector<double>     vA = Notation::Kelvin::to_vector(A);
  const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);
  Vector<double>           vB(mC.m());
  mC.vmult(vB, vA);

  typedef typename std::decay<decltype(A)>::type InpVecType;
  typedef typename std::decay<decltype(B)>::type ResVecType;
  typedef typename std::decay<decltype(C)>::type InpMatType;
  const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
  const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
  const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

  std::cout << "Scalar" << std::endl;
  std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
  std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
  std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

  Assert(std::abs(A_conv - A) < 1e-12,
         ExcMessage("Different result for input A"));
  Assert(std::abs(C_conv - C) < 1e-12,
         ExcMessage("Different result for input C"));
  Assert(std::abs(B_conv - B) < 1e-12,
         ExcMessage("Different result for output B"));

  Assert(std::abs(vA.l2_norm() - std::abs(A)) < 1e-12,
         ExcMessage("Different norm for input A"));
  Assert(std::abs(mC.frobenius_norm() - std::abs(C)) < 1e-12,
         ExcMessage("Different norm for input C"));
  Assert(std::abs(vB.l2_norm() - std::abs(B)) < 1e-12,
         ExcMessage("Different norm for output B"));
}

template <int dim>
void
test_rank_0_tensors()
{
  const Tensor<0, dim, double> A = 5;
  const Tensor<0, dim, double> C = 12;
  const Tensor<0, dim, double> B = C * A;

  const Vector<double>     vA = Notation::Kelvin::to_vector(A);
  const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);
  Vector<double>           vB(mC.m());
  mC.vmult(vB, vA);

  typedef typename std::decay<decltype(A)>::type InpVecType;
  typedef typename std::decay<decltype(B)>::type ResVecType;
  typedef typename std::decay<decltype(C)>::type InpMatType;
  const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
  const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
  const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

  std::cout << "Rank 0" << std::endl;
  std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
  std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
  std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

  Assert((A_conv - A).norm() < 1e-12,
         ExcMessage("Different result for input A"));
  Assert((C_conv - C).norm() < 1e-12,
         ExcMessage("Different result for input C"));
  Assert((B_conv - B).norm() < 1e-12,
         ExcMessage("Different result for output B"));

  Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
         ExcMessage("Different norm for input A"));
  Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
         ExcMessage("Different norm for input C"));
  Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
         ExcMessage("Different norm for output B"));
}

template <int dim>
void
test_rank_1_2_tensors()
{
  // Non-symmetric tensor
  {
    Tensor<1, dim, double> A;
    Tensor<2, dim, double> C;
    initialize(A);
    initialize(C);
    const Tensor<1, dim, double> B = C * A;

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);
    Vector<double>           vB(mC.m());
    mC.vmult(vB, vA);

    typedef typename std::decay<decltype(A)>::type InpVecType;
    typedef typename std::decay<decltype(B)>::type ResVecType;
    typedef typename std::decay<decltype(C)>::type InpMatType;
    const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
    const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
    const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

    std::cout << "Rank 1 (non-symm)" << std::endl;
    std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
    std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
    std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

    Assert((A_conv - A).norm() < 1e-12,
           ExcMessage("Different result for input A"));
    Assert((C_conv - C).norm() < 1e-12,
           ExcMessage("Different result for input C"));
    Assert((B_conv - B).norm() < 1e-12,
           ExcMessage("Different result for output B"));

    Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
           ExcMessage("Different norm for input A"));
    Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
           ExcMessage("Different norm for input C"));
    Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
           ExcMessage("Different norm for output B"));
  }

  // Symmetric tensor
  {
    Tensor<1, dim, double>          A;
    SymmetricTensor<2, dim, double> C;
    initialize(A);
    initialize(C);
    const Tensor<1, dim, double> B = C * A;

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);
    Vector<double>           vB(mC.m());
    mC.vmult(vB, vA);

    typedef typename std::decay<decltype(A)>::type InpVecType;
    typedef typename std::decay<decltype(B)>::type ResVecType;
    typedef typename std::decay<decltype(C)>::type InpMatType;
    const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
    const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
    const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

    std::cout << "Rank 1 (symm)" << std::endl;
    std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
    std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
    std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

    Assert((A_conv - A).norm() < 1e-12,
           ExcMessage("Different result for input A"));
    Assert((C_conv - C).norm() < 1e-12,
           ExcMessage("Different result for input C"));
    Assert((B_conv - B).norm() < 1e-12,
           ExcMessage("Different result for output B"));

    Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
           ExcMessage("Different norm for input A"));
    Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
           ExcMessage("Different norm for input C"));
    Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
           ExcMessage("Different norm for output B"));
  }
}

template <int dim>
void
test_rank_2_4_tensors()
{
  // Non-symmetric tensor
  {
    Tensor<2, dim, double> A;
    Tensor<4, dim, double> C;
    initialize(A);
    initialize(C);
    const Tensor<2, dim, double> B = double_contract<2, 0, 3, 1>(C, A);

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);
    Vector<double>           vB(mC.m());
    mC.vmult(vB, vA);

    typedef typename std::decay<decltype(A)>::type InpVecType;
    typedef typename std::decay<decltype(B)>::type ResVecType;
    typedef typename std::decay<decltype(C)>::type InpMatType;
    const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
    const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
    const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

    std::cout << "Rank 2 (non-symm)" << std::endl;
    std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
    std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
    std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

    Assert((A_conv - A).norm() < 1e-12,
           ExcMessage("Different result for input A"));
    Assert((C_conv - C).norm() < 1e-12,
           ExcMessage("Different result for input C"));
    Assert((B_conv - B).norm() < 1e-12,
           ExcMessage("Different result for output B"));

    Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
           ExcMessage("Different norm for input A"));
    Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
           ExcMessage("Different norm for input C"));
    Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
           ExcMessage("Different norm for output B"));
  }

  // Symmetric tensor
  {
    SymmetricTensor<2, dim, double> A;
    SymmetricTensor<4, dim, double> C;
    initialize(A);
    initialize(C);
    const SymmetricTensor<2, dim, double> B = C * A;

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);
    Vector<double>           vB(mC.m());
    mC.vmult(vB, vA);

    typedef typename std::decay<decltype(A)>::type InpVecType;
    typedef typename std::decay<decltype(B)>::type ResVecType;
    typedef typename std::decay<decltype(C)>::type InpMatType;
    const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
    const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
    const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

    std::cout << "Rank 2 (symm)" << std::endl;
    std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
    std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
    std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

    Assert((A_conv - A).norm() < 1e-12,
           ExcMessage("Different result for input A"));
    Assert((C_conv - C).norm() < 1e-12,
           ExcMessage("Different result for input C"));
    Assert((B_conv - B).norm() < 1e-12,
           ExcMessage("Different result for output B"));

    Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
           ExcMessage("Different norm for input A"));
    Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
           ExcMessage("Different norm for input C"));
    Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
           ExcMessage("Different norm for output B"));
  }

  // Non-symmetric tensor from symmetric tensor
  {
    SymmetricTensor<2, dim, double> A_symm;
    SymmetricTensor<4, dim, double> C_symm;
    initialize(A_symm);
    initialize(C_symm);

    Tensor<2, dim, double>       A(A_symm);
    Tensor<4, dim, double>       C(C_symm);
    const Tensor<2, dim, double> B = double_contract<2, 0, 3, 1>(C, A);

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);
    Vector<double>           vB(mC.m());
    mC.vmult(vB, vA);

    typedef typename std::decay<decltype(A)>::type InpVecType;
    typedef typename std::decay<decltype(B)>::type ResVecType;
    typedef typename std::decay<decltype(C)>::type InpMatType;
    const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
    const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
    const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

    std::cout << "Rank 2 (non-symm from symm)" << std::endl;
    std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
    std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
    std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

    Assert((A_conv - A).norm() < 1e-12,
           ExcMessage("Different result for input A"));
    Assert((C_conv - C).norm() < 1e-12,
           ExcMessage("Different result for input C"));
    Assert((B_conv - B).norm() < 1e-12,
           ExcMessage("Different result for output B"));

    Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
           ExcMessage("Different norm for input A"));
    Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
           ExcMessage("Different norm for input C"));
    Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
           ExcMessage("Different norm for output B"));
  }
}

template <int dim>
void
test_rank_3_tensors()
{
  // Non-symmetric tensor: Version 1
  {
    Tensor<2, dim, double> A;
    Tensor<3, dim, double> C;
    initialize(A);
    initialize(C);
    const Tensor<1, dim, double> B =
      double_contract<0, 0, 1, 1>(C,
                                  A); // This implies that a Tvmult is necessary

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);
    Vector<double>           vB(mC.n()); // Note result size
    mC.Tvmult(vB, vA);                   // Note transpose vmult

    typedef typename std::decay<decltype(A)>::type InpVecType;
    typedef typename std::decay<decltype(B)>::type ResVecType;
    typedef typename std::decay<decltype(C)>::type InpMatType;
    const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
    const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
    const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

    std::cout << "Rank 3 (non-symm 1)" << std::endl;
    std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
    std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
    std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

    Assert((A_conv - A).norm() < 1e-12,
           ExcMessage("Different result for input A"));
    Assert((C_conv - C).norm() < 1e-12,
           ExcMessage("Different result for input C"));
    Assert((B_conv - B).norm() < 1e-12,
           ExcMessage("Different result for output B"));

    Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
           ExcMessage("Different norm for input A"));
    Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
           ExcMessage("Different norm for input C"));
    Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
           ExcMessage("Different norm for output B"));
  }

  // Non-symmetric tensor: Version 2
  {
    Tensor<2, dim, double> A;
    Tensor<3, dim, double> C;
    initialize(A);
    initialize(C);
    const Tensor<1, dim, double> B = double_contract<1, 0, 2, 1>(
      C, A); // This implies that a standard vmult is necessary

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mC =
      Notation::Kelvin::to_matrix<dim, Tensor<1, dim>, Tensor<2, dim>>(
        C); // Define subtensor representation Tensor<1,dim> \otimes
            // Tensor<2,dim>
    Vector<double> vB(mC.m()); // Note result size
    mC.vmult(vB, vA);          // Note transpose vmult

    typedef typename std::decay<decltype(A)>::type InpVecType;
    typedef typename std::decay<decltype(B)>::type ResVecType;
    typedef typename std::decay<decltype(C)>::type InpMatType;
    const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
    const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
    const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

    std::cout << "Rank 3 (non-symm 2)" << std::endl;
    std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
    std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
    std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

    Assert((A_conv - A).norm() < 1e-12,
           ExcMessage("Different result for input A"));
    Assert((C_conv - C).norm() < 1e-12,
           ExcMessage("Different result for input C"));
    Assert((B_conv - B).norm() < 1e-12,
           ExcMessage("Different result for output B"));

    Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
           ExcMessage("Different norm for input A"));
    Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
           ExcMessage("Different norm for input C"));
    Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
           ExcMessage("Different norm for output B"));
  }

  // Symmetric tensor: Version 1
  {
    SymmetricTensor<2, dim, double> A;
    Tensor<3, dim, double>          C;
    initialize(A);
    initialize(C, true); // Specialised constructor
    const Tensor<2, dim, double> A_ns(A);
    const Tensor<1, dim, double> B = double_contract<0, 0, 1, 1>(
      C, A_ns); // This implies that a Tvmult is necessary

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mC =
      Notation::Kelvin::to_matrix<dim, SymmetricTensor<2, dim>, Tensor<1, dim>>(
        C); // Define subtensor representation SymmetricTensor<2,dim> \otimes
            // Tensor<1,dim>
    Vector<double> vB(mC.n()); // Note result size
    mC.Tvmult(vB, vA);         // Note transpose vmult

    typedef typename std::decay<decltype(A)>::type InpVecType;
    typedef typename std::decay<decltype(B)>::type ResVecType;
    typedef typename std::decay<decltype(C)>::type InpMatType;
    const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
    const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
    const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

    std::cout << "Rank 3 (symm 1)" << std::endl;
    std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
    std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
    std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

    Assert((A_conv - A).norm() < 1e-12,
           ExcMessage("Different result for input A"));
    Assert((C_conv - C).norm() < 1e-12,
           ExcMessage("Different result for input C"));
    Assert((B_conv - B).norm() < 1e-12,
           ExcMessage("Different result for output B"));

    Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
           ExcMessage("Different norm for input A"));
    Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
           ExcMessage("Different norm for input C"));
    Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
           ExcMessage("Different norm for output B"));
  }

  // Symmetric tensor: Version 2
  {
    SymmetricTensor<2, dim, double> A;
    Tensor<3, dim, double>          C;
    initialize(A);
    initialize(C, false); // Specialised constructor
    const Tensor<2, dim, double> A_ns(A);
    const Tensor<1, dim, double> B = double_contract<1, 0, 2, 1>(
      C, A_ns); // This implies that a standard vmult is necessary

    const Vector<double>     vA = Notation::Kelvin::to_vector(A);
    const FullMatrix<double> mC =
      Notation::Kelvin::to_matrix<dim, Tensor<1, dim>, SymmetricTensor<2, dim>>(
        C); // Define subtensor representation Tensor<1,dim> \otimes
            // SymmetricTensor<2,dim>
    Vector<double> vB(mC.m()); // Note result size
    mC.vmult(vB, vA);          // Note transpose vmult

    typedef typename std::decay<decltype(A)>::type InpVecType;
    typedef typename std::decay<decltype(B)>::type ResVecType;
    typedef typename std::decay<decltype(C)>::type InpMatType;
    const auto A_conv = Notation::Kelvin::to_tensor<InpVecType>(vA);
    const auto B_conv = Notation::Kelvin::to_tensor<ResVecType>(vB);
    const auto C_conv = Notation::Kelvin::to_tensor<InpMatType>(mC);

    std::cout << "Rank 3 (symm 2)" << std::endl;
    std::cout << "A: " << A << "  A_conv: " << A_conv << std::endl;
    std::cout << "B: " << B << "  B_conv: " << B_conv << std::endl;
    std::cout << "C: " << C << "  C_conv: " << C_conv << std::endl;

    Assert((A_conv - A).norm() < 1e-12,
           ExcMessage("Different result for input A"));
    Assert((C_conv - C).norm() < 1e-12,
           ExcMessage("Different result for input C"));
    Assert((B_conv - B).norm() < 1e-12,
           ExcMessage("Different result for output B"));

    Assert(std::abs(vA.l2_norm() - A.norm()) < 1e-12,
           ExcMessage("Different norm for input A"));
    Assert(std::abs(mC.frobenius_norm() - C.norm()) < 1e-12,
           ExcMessage("Different norm for input C"));
    Assert(std::abs(vB.l2_norm() - B.norm()) < 1e-12,
           ExcMessage("Different norm for output B"));
  }
}

template <int dim>
void
test_tensors()
{
  test_scalars<dim>();
  test_rank_0_tensors<dim>();
  test_rank_1_2_tensors<dim>();
  test_rank_2_4_tensors<dim>();
  test_rank_3_tensors<dim>();
}

int
main()
{
  initlog();

  test_tensors<2>();
  test_tensors<3>();

  deallog << "OK" << std::endl;
}
