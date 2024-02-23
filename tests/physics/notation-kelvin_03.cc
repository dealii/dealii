// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check the scaling factors used during conversion between
// Tensor and Kelvin notation

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/physics/notation.h>

#include "../tests.h"


using namespace dealii::Physics;

template <int dim, typename Number>
void
initialize(Tensor<0, dim, Number> &x)
{
  x = 1.0;
}

template <int dim, typename Number>
void
initialize(Tensor<1, dim, Number> &x)
{
  for (unsigned int i = 0; i < x.n_independent_components; ++i)
    x[i] = 1.0;
}

template <int dim, typename Number>
void
initialize(Tensor<2, dim, Number> &x)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        x[i][j] = 1.0;
      }
}

template <int dim, typename Number>
void
initialize(SymmetricTensor<2, dim, Number> &x)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      {
        x[i][j] = 1.0;
      }
}

template <int dim, typename Number>
void
initialize(Tensor<3, dim, Number> &x)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        {
          x[i][j][k] = 1.0;
        }
}

template <int dim, typename Number>
void
initialize(Tensor<4, dim, Number> &x)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          {
            x[i][j][k][l] = 1.0;
          }
}

template <int dim, typename Number>
void
initialize(SymmetricTensor<4, dim, Number> &x)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          {
            x[i][j][k][l] = 1.0;
          }
}

template <int dim>
void
test_scalars()
{
  const double A = 1.0;
  const double C = 1.0;

  const Vector<double>     vA = Notation::Kelvin::to_vector(A);
  const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);

  deallog.push("Scalar");
  {
    deallog << "Vector" << std::endl;
    vA.print(deallog.get_file_stream());
    deallog << "Matrix" << std::endl;
    mC.print_formatted(deallog.get_file_stream());
  }
  deallog.pop();
}

template <int rank, int dim>
void
test_tensors_012()
{
  static_assert((rank == 0 || rank == 1 || rank == 2),
                "Must be rank 0, 1 or 2");

  Tensor<rank, dim, double> A;
  Tensor<rank, dim, double> C;
  initialize(A);
  initialize(C);

  const Vector<double>     vA = Notation::Kelvin::to_vector(A);
  const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);

  const std::string name = std::string("Rank-") +
                           Utilities::int_to_string(rank) +
                           std::string(" tensor");
  deallog.push(name);
  {
    deallog << "Vector" << std::endl;
    vA.print(deallog.get_file_stream());
    deallog << "Matrix" << std::endl;
    mC.print_formatted(deallog.get_file_stream());
  }
  deallog.pop();
}

template <int dim>
void
test_tensor_3()
{
  // Non-symmetric tensor: Version 1
  {
    Tensor<3, dim, double> C;
    initialize(C);

    const FullMatrix<double> mC =
      Notation::Kelvin::to_matrix<dim, Tensor<1, dim>, Tensor<2, dim>>(C);

    deallog.push("Rank-3 tensor: NSymm 1");
    {
      deallog << "Matrix" << std::endl;
      mC.print_formatted(deallog.get_file_stream());
    }
    deallog.pop();
  }

  // Non-symmetric tensor: Version 2
  {
    Tensor<3, dim, double> C;
    initialize(C);

    const FullMatrix<double> mC =
      Notation::Kelvin::to_matrix<dim, Tensor<2, dim>, Tensor<1, dim>>(C);

    deallog.push("Rank-3 tensor: NSymm 2");
    {
      deallog << "Matrix" << std::endl;
      mC.print_formatted(deallog.get_file_stream());
    }
    deallog.pop();
  }

  // Symmetric tensor: Version 1
  {
    Tensor<3, dim, double> C;
    initialize(C);

    const FullMatrix<double> mC =
      Notation::Kelvin::to_matrix<dim, Tensor<1, dim>, SymmetricTensor<2, dim>>(
        C);

    deallog.push("Rank-3 tensor: Symm 1");
    {
      deallog << "Matrix" << std::endl;
      mC.print_formatted(deallog.get_file_stream());
    }
    deallog.pop();
  }

  // Symmetric tensor: Version 2
  {
    Tensor<3, dim, double> C;
    initialize(C);

    const FullMatrix<double> mC =
      Notation::Kelvin::to_matrix<dim, SymmetricTensor<2, dim>, Tensor<1, dim>>(
        C);

    deallog.push("Rank-3 tensor: Symm 2");
    {
      deallog << "Matrix" << std::endl;
      mC.print_formatted(deallog.get_file_stream());
    }
    deallog.pop();
  }
}

template <int dim>
void
test_tensor_4()
{
  Tensor<4, dim, double> C;
  initialize(C);

  const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);

  deallog.push("Rank-4 tensor");
  {
    deallog << "Matrix" << std::endl;
    mC.print_formatted(deallog.get_file_stream());
  }
  deallog.pop();
}

template <int dim>
void
test_symmetric_tensor_2()
{
  SymmetricTensor<2, dim, double> A;
  SymmetricTensor<2, dim, double> C;
  initialize(A);
  initialize(C);

  const Vector<double>     vA = Notation::Kelvin::to_vector(A);
  const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);

  deallog.push("Rank-2 symm tensor");
  {
    deallog << "Vector" << std::endl;
    vA.print(deallog.get_file_stream());
    deallog << "Matrix" << std::endl;
    mC.print_formatted(deallog.get_file_stream());
  }
  deallog.pop();
}

template <int dim>
void
test_symmetric_tensor_4()
{
  SymmetricTensor<4, dim, double> C;
  initialize(C);

  const FullMatrix<double> mC = Notation::Kelvin::to_matrix(C);

  deallog.push("Rank-4 symm tensor");
  {
    deallog << "Matrix" << std::endl;
    mC.print_formatted(deallog.get_file_stream());
  }
  deallog.pop();
}

template <int dim>
void
test_tensors()
{
  deallog.push(Utilities::int_to_string(dim));
  {
    test_scalars<dim>();
    test_tensors_012<0, dim>();
    test_tensors_012<1, dim>();
    test_tensors_012<2, dim>();
    test_tensor_3<dim>();
    test_tensor_4<dim>();
    test_symmetric_tensor_2<dim>();
    test_symmetric_tensor_4<dim>();

    deallog << "OK" << std::endl;
  }
  deallog.pop();
}

int
main()
{
  initlog();

  test_tensors<2>();
  test_tensors<3>();

  deallog << "OK" << std::endl;
}
