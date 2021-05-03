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


// Check the contract3 function involving SymmetricTensors

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include "../tests.h"


// Although not strictly required, in this test it is expected that T1 or T3
// (or both) are a symmetric tensor
template <int rank1,
          int rank2,
          int rank3,
          int dim,
          typename number,
          template <int, int, typename> class T1,
          template <int, int, typename> class T3>
void
test_symm_tensor_contract_3(const T1<rank1, dim, number> &    l,
                            const Tensor<rank2, dim, number> &m,
                            const T3<rank3, dim, number> &    r)
{
  const double res1 = contract3(l, m, r);
  const double res2 = contract3(static_cast<Tensor<rank1, dim>>(l),
                                m,
                                static_cast<Tensor<rank3, dim>>(r));
  Assert(std::abs(res1 - res2) < 1e-12,
         ExcMessage("Result from symmetric tensor contract3 is incorrect."));
}

int
main()
{
  initlog();
  deallog << std::setprecision(5);

  const int               dim = 3;
  Tensor<1, dim>          v1;
  Tensor<2, dim>          t1;
  SymmetricTensor<2, dim> s1, s2;
  Tensor<3, dim>          T1;
  Tensor<4, dim>          H1;

  for (unsigned int i = 0; i < dim; ++i)
    {
      v1[i] = 1 + i;
      for (unsigned int j = 0; j < dim; ++j)
        {
          t1[i][j] = 1 + j + i * dim;
          s1[i][j] = 1 + j + i * dim;
          s2[i][j] = 2 + (j + i) * dim;
          for (unsigned int k = 0; k < dim; ++k)
            {
              T1[i][j][k] = 1 + k + j * dim + i * dim * dim;
              for (unsigned int l = 0; l < dim; ++l)
                H1[i][j][k][l] =
                  1 + l + k * dim + j * dim * dim + i * dim * dim * dim;
            }
        }
    }

  // Vector + SymmetricTensor
  test_symm_tensor_contract_3(v1, T1, s1);
  test_symm_tensor_contract_3(s1, T1, v1);

  // Tensor + SymmetricTensor
  test_symm_tensor_contract_3(t1, H1, s1);
  test_symm_tensor_contract_3(s1, H1, t1);

  // SymmetricTensor + SymmetricTensor
  test_symm_tensor_contract_3(s1, H1, s2);

  deallog << "All OK" << std::endl;
}
