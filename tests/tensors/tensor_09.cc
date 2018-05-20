// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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

// check generic contract variants tha allow to specify the indices which
// will be contracted.

#include "../tests.h"
#include <deal.II/base/tensor.h>

int
main()
{
  initlog();

  Tensor<1, 3, int> rank1;
  rank1[0] = 1;
  rank1[1] = 2;
  rank1[2] = 4;

  Tensor<2, 3, int> rank2;
  rank2[0] = 4 * rank1;
  rank2[1] = 2 * rank1;
  rank2[2] = 1 * rank1;

  Tensor<3, 3, int> rank3;
  rank3[0] = 7 * rank2;
  rank3[1] = 6 * rank2;
  rank3[2] = 9 * rank2;

  Tensor<4, 3, int> rank4;
  rank4[0] = 4 * rank3;
  rank4[1] = 3 * rank3;
  rank4[2] = 2 * rank3;

  // rank3 * rank1 over all possible indices:
  {
    Tensor<2, 3, int> result;

    deallog << contract<0, 0>(rank3, rank1) << std::endl;
    contract(result, rank3, 1, rank1);
    deallog << result << std::endl << std::endl;

    deallog << contract<1, 0>(rank3, rank1) << std::endl;
    contract(result, rank3, 2, rank1);
    deallog << result << std::endl << std::endl;

    deallog << contract<2, 0>(rank3, rank1) << std::endl;
    contract(result, rank3, 3, rank1);
    deallog << result << std::endl << std::endl;
  }

  // rank2 * rank2 over all possible indices:
  {
    Tensor<2, 3, int> result;

    deallog << contract<0, 0>(rank2, rank2) << std::endl;
    contract(result, rank2, 1, rank2, 1);
    deallog << result << std::endl << std::endl;

    deallog << contract<1, 0>(rank2, rank2) << std::endl;
    contract(result, rank2, 2, rank2, 1);
    deallog << result << std::endl << std::endl;

    deallog << contract<0, 1>(rank2, rank2) << std::endl;
    contract(result, rank2, 1, rank2, 2);
    deallog << result << std::endl << std::endl;

    deallog << contract<1, 1>(rank2, rank2) << std::endl;
    contract(result, rank2, 2, rank2, 2);
    deallog << result << std::endl << std::endl;
  }

  // rank3 * rank2 over all possible indices:
  {
    Tensor<3, 3, int> result;

    deallog << contract<0, 0>(rank3, rank2) << std::endl;
    contract(result, rank3, 1, rank2, 1);
    deallog << result << std::endl << std::endl;

    deallog << contract<1, 0>(rank3, rank2) << std::endl;
    contract(result, rank3, 2, rank2, 1);
    deallog << result << std::endl << std::endl;

    deallog << contract<2, 0>(rank3, rank2) << std::endl;
    contract(result, rank3, 3, rank2, 1);
    deallog << result << std::endl << std::endl;

    deallog << contract<0, 1>(rank3, rank2) << std::endl;
    contract(result, rank3, 1, rank2, 2);
    deallog << result << std::endl << std::endl;

    deallog << contract<1, 1>(rank3, rank2) << std::endl;
    contract(result, rank3, 2, rank2, 2);
    deallog << result << std::endl << std::endl;

    deallog << contract<2, 1>(rank3, rank2) << std::endl;
    contract(result, rank3, 3, rank2, 2);
    deallog << result << std::endl << std::endl;
  }

  // rank4 ** rank2, double contraction:
  {
    Tensor<2, 3, int> result;

    deallog << double_contract<2, 0, 3, 1>(rank4, rank2) << std::endl;
    deallog << double_contract<3, 1, 2, 0>(rank4, rank2) << std::endl;
    deallog << double_contract<0, 2, 1, 3>(rank2, rank4) << std::endl;
    deallog << double_contract<1, 3, 0, 2>(rank2, rank4) << std::endl;

    double_contract(result, rank4, rank2);
    deallog << result << std::endl << std::endl;
  }
}
