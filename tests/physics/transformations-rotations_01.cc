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


// test rotation matrix definitions

#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/physics/transformations.h>

#include "../tests.h"


using namespace dealii::Physics;


void
test_rotation_matrix_3d_z_axis(const double angle)
{
  Tensor<2, 3>       R_z = unit_symmetric_tensor<3>();
  const Tensor<2, 2> R_2d =
    Transformations::Rotations::rotation_matrix_2d(angle);
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      R_z[i][j] = R_2d[i][j];

  Assert(std::abs(determinant(R_z) - 1.0) < 1e-9,
         ExcMessage("Rodrigues rotation matrix determinant is not unity"));
  const Tensor<2, 3> R =
    Transformations::Rotations::rotation_matrix_3d(Point<3>({0, 0, 1}), angle);
  Assert(std::abs(determinant(R) - 1.0) < 1e-9,
         ExcMessage("Rotation matrix determinant is not unity"));

  Assert((transpose(R) * R - unit_symmetric_tensor<3>()).norm() < 1e-9,
         ExcMessage("Matrix is not a rotation matrix"));
  Assert((R - R_z).norm() < 1e-12,
         ExcMessage("Incorrect computation of R in 3d"));
}

void
test_rotation_matrix_3d(const Point<3> &axis, const double angle)
{
  // http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
  // http://en.wikipedia.org/wiki/Rotation_matrix
  // NOTE: Angle in radians
  const Tensor<1, 3> u        = axis / axis.norm(); // Ensure unit vector
  const Tensor<2, 3> u_dyad_u = outer_product(u, u);
  const double       u_skew_array[3][3] = {{0.0, -u[2], u[1]},
                                     {u[2], 0.0, -u[0]},
                                     {-u[1], u[0], 0.0}};

  const Tensor<2, 3> R_rodrigues =
    u_dyad_u +
    std::cos(angle) *
      (static_cast<Tensor<2, 3>>(unit_symmetric_tensor<3>()) - u_dyad_u) +
    std::sin(angle) * Tensor<2, 3>(u_skew_array);


  Assert(std::abs(determinant(R_rodrigues) - 1.0) < 1e-9,
         ExcMessage("Rodrigues rotation matrix determinant is not unity"));
  const Tensor<2, 3> R =
    Transformations::Rotations::rotation_matrix_3d(axis, angle);
  Assert(std::abs(determinant(R) - 1.0) < 1e-9,
         ExcMessage("Rotation matrix determinant is not unity"));

  Assert((transpose(R) * R - unit_symmetric_tensor<3>()).norm() < 1e-9,
         ExcMessage("Matrix is not a rotation matrix"));
  Assert((R - R_rodrigues).norm() < 1e-12,
         ExcMessage("Incorrect computation of R in 3d"));
}

Point<3>
normalise(const Point<3> &p)
{
  Assert(p.norm() > 0.0, ExcMessage("Point vector has zero norm"));
  return p / p.norm();
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  const double deg_to_rad = numbers::PI / 180.0;

  // 2-d
  {
    const double       angle = 90.0 * deg_to_rad;
    const Tensor<1, 2> in({1, 0});
    const Tensor<2, 2> R =
      Transformations::Rotations::rotation_matrix_2d(angle);
    Assert((transpose(R) * R - unit_symmetric_tensor<2>()).norm() < 1e-9,
           ExcMessage("Matrix is not a rotation matrix"));
    Assert(std::abs(determinant(R) - 1.0) < 1e-9,
           ExcMessage("Rotation matrix determinant is not unity"));
    const Tensor<1, 2> out = R * in;
    Assert((out - Tensor<1, 2>({0, 1})).norm() < 1e-12,
           ExcMessage("Incorrect computation of 90 degree R in 2d"));
  }
  {
    const double       angle = 135.0 * deg_to_rad;
    const Tensor<1, 2> in({1, 0});
    const Tensor<2, 2> R =
      Transformations::Rotations::rotation_matrix_2d(angle);
    Assert((transpose(R) * R - unit_symmetric_tensor<2>()).norm() < 1e-9,
           ExcMessage("Matrix is not a rotation matrix"));
    Assert(std::abs(determinant(R) - 1.0) < 1e-9,
           ExcMessage("Rotation matrix determinant is not unity"));
    const Tensor<1, 2> out = R * in;
    Assert((out - Tensor<1, 2>({-1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0)}))
               .norm() < 1e-12,
           ExcMessage("Incorrect computation of 135 degree R in 2d"));
  }
  {
    const double       angle = 240.0 * deg_to_rad;
    const Tensor<1, 2> in({1, 0});
    const Tensor<2, 2> R =
      Transformations::Rotations::rotation_matrix_2d(angle);
    Assert((transpose(R) * R - unit_symmetric_tensor<2>()).norm() < 1e-9,
           ExcMessage("Matrix is not a rotation matrix"));
    Assert(std::abs(determinant(R) - 1.0) < 1e-9,
           ExcMessage("Rotation matrix determinant is not unity"));
    const Tensor<1, 2> out = R * in;
    Assert((out - Tensor<1, 2>({-0.5, -std::sqrt(3.0) / 2.0})).norm() < 1e-12,
           ExcMessage("Incorrect computation of 240 degree R in 2d"));
  }

  // 3-d
  test_rotation_matrix_3d_z_axis(90.0 * deg_to_rad);
  test_rotation_matrix_3d_z_axis(45.0 * deg_to_rad);
  test_rotation_matrix_3d_z_axis(60.0 * deg_to_rad);

  test_rotation_matrix_3d(normalise(Point<3>({1, 1, 1})), 90.0 * deg_to_rad);
  test_rotation_matrix_3d(normalise(Point<3>({0, 2, 1})), 45.0 * deg_to_rad);
  test_rotation_matrix_3d(normalise(Point<3>({-1, 3, 2})), 60.0 * deg_to_rad);

  deallog << "OK" << std::endl;
}
