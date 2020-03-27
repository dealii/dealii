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

// Test vectorization capabilities of the
// dealii::Physics::Elasticity::kinematics quantities.

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>

#include "../tests.h"


using namespace dealii::Physics;
using namespace dealii::Physics::Elasticity;

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  const int dim = 3;

  Tensor<2, dim, VectorizedArray<double>> grad_u;

  // grad_u = [2 -1 0; -1 2 -1; 0 -1 2]
  // which is a random s.p.d. matrix -> F will be s.p.d
  // -> det(F) > 0 which has to be satisfied for
  // F to make any sense (no negative volume).
  grad_u[0][0] = 2.0;
  grad_u[0][1] = -1.0;
  grad_u[0][2] = 0.0;
  grad_u[1][0] = -1.0;
  grad_u[1][1] = 2.0;
  grad_u[1][2] = -1.0;
  grad_u[2][0] = 0.0;
  grad_u[2][1] = -1.0;
  grad_u[2][2] = 2.0;

  // Scale the gradients along the vectorization-index so that each grad_u[v] is
  // unique.
  for (unsigned int v = 0; v < VectorizedArray<double>::size(); v++)
    for (unsigned int i = 0; i < dim; i++)
      for (unsigned int j = 0; j < dim; j++)
        grad_u[i][j][v] *= (v + 1);

  Tensor<2, dim, VectorizedArray<double>> F_solution;
  F_solution = grad_u;
  for (unsigned int i = 0; i < dim; i++)
    F_solution[i][i] = F_solution[i][i] + 1.0;

  Tensor<2, dim, VectorizedArray<double>> F_test = Kinematics::F(grad_u);

  // You can't use .norm() on some difference-tensor of the two so we compare
  // element-wise!
  for (unsigned int v = 0; v < VectorizedArray<double>::size(); v++)
    for (unsigned int i = 0; i < dim; i++)
      for (unsigned int j = 0; j < dim; j++)
        if (F_solution[i][j][v] - F_test[i][j][v] != 0.0)
          deallog << "Not OK" << std::endl;

  SymmetricTensor<2, dim, VectorizedArray<double>> E_solution;
  E_solution =
    0.5 * (symmetrize(transpose(F_solution) * F_solution) -
           static_cast<SymmetricTensor<2, dim, VectorizedArray<double>>>(
             StandardTensors<dim>::I));

  SymmetricTensor<2, dim, VectorizedArray<double>> E_test =
    Kinematics::E(F_test);

  for (unsigned int v = 0; v < VectorizedArray<double>::size(); v++)
    for (unsigned int i = 0; i < dim; i++)
      for (unsigned int j = 0; j < dim; j++)
        if (E_test[i][j][v] - E_solution[i][j][v] != 0.0)
          deallog << "Not OK" << std::endl;

  Tensor<2, dim, VectorizedArray<double>> F_iso_solution;
  F_iso_solution = std::pow(determinant(F_solution), -1.0 / dim) * F_solution;

  Tensor<2, dim, VectorizedArray<double>> F_iso_test;
  F_iso_test = Kinematics::F_iso(F_test);

  for (unsigned int v = 0; v < VectorizedArray<double>::size(); v++)
    for (unsigned int i = 0; i < dim; i++)
      for (unsigned int j = 0; j < dim; j++)
        if (F_iso_test[i][j][v] - F_iso_solution[i][j][v] != 0.0)
          deallog << "Not OK" << std::endl;

  deallog << "OK" << std::endl;
}
