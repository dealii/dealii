// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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


// Check the inverse of a rank-2 Tensor

// Equivalent matlab script:
/*
printf ("Symmetric Tensor dim 1\n")
t0 = [2]
inv(t0)

printf ("Symmetric Tensor dim 2\n")
t1 = [2 -1; -1 1.5]
inv(t1)

printf ("Symmetric Tensor dim 3\n")
t2 = [2   1    1.5;
      1   1.5  0.25;
      1.5 0.25 1.25]
inv(t2)
*/

// Symmetric Tensor dim 1
// t0 =  2
// ans =  0.50000
// Symmetric Tensor dim 2
// t1 =
//
//    2.0000  -1.0000
//   -1.0000   1.5000
//
// ans =
//
//    0.75000   0.50000
//    0.50000   1.00000
//
// Symmetric Tensor dim 3
// t2 =
//
//    2.00000   1.00000   1.50000
//    1.00000   1.50000   0.25000
//    1.50000   0.25000   1.25000
//
// ans =
//
//   -7.2500   3.5000   8.0000
//    3.5000  -1.0000  -4.0000
//    8.0000  -4.0000  -8.0000

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(5);

  deallog << "Symmetric Tensor dim 1" << std::endl;
  SymmetricTensor<2, 1> t1;
  t1[0][0] = 2.0;
  deallog << invert(t1) << std::endl;
  Assert((static_cast<Tensor<2, 1>>(invert(t1)) *
            static_cast<Tensor<2, 1>>(t1) -
          unit_symmetric_tensor<1>())
             .norm() < 1e-12,
         ExcMessage("Dim 1 inverse symmetric tensor definition is incorrect"));

  deallog << "Symmetric Tensor dim 2" << std::endl;
  SymmetricTensor<2, 2> t2;
  t2[0][0] = 2.0;
  t2[0][1] = 1.0;
  t2[1][1] = 1.5;
  deallog << invert(t2) << std::endl;
  Assert((static_cast<Tensor<2, 2>>(invert(t2)) *
            static_cast<Tensor<2, 2>>(t2) -
          unit_symmetric_tensor<2>())
             .norm() < 1e-12,
         ExcMessage("Dim 2 inverse symmetric tensor definition is incorrect"));

  deallog << "Symmetric Tensor dim 3" << std::endl;
  SymmetricTensor<2, 3> t3;
  t3[0][0] = 2.0;
  t3[0][1] = 1.0;
  t3[0][2] = 1.5;
  t3[1][1] = 1.5;
  t3[1][2] = 0.25;
  t3[2][2] = 1.25;
  deallog << invert(t3) << std::endl;
  Assert((static_cast<Tensor<2, 3>>(invert(t3)) *
            static_cast<Tensor<2, 3>>(t3) -
          unit_symmetric_tensor<3>())
             .norm() < 1e-12,
         ExcMessage("Dim 3 inverse symmetric tensor definition is incorrect"));

  deallog << "OK" << std::endl;
}
