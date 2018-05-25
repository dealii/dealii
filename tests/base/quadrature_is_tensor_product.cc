// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// check the boolean is_tensor_product for all the quadrature classes


#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

template <int dim>
void
print_is_tensor_product()
{
  Quadrature<dim> q_0(1);
  deallog << "Quadrature<" << dim << ">: " << q_0.is_tensor_product()
          << std::endl;

  QIterated<dim> q_1(QGauss<1>(1), 1);
  deallog << "QIterated<" << dim << ">: " << q_1.is_tensor_product()
          << std::endl;

  QGauss<dim> q_2(1);
  deallog << "QGauss<" << dim << ">: " << q_2.is_tensor_product() << std::endl;

  QGaussLobatto<dim> q_3(2);
  deallog << "QGaussLobatto<" << dim << ">: " << q_3.is_tensor_product()
          << std::endl;

  QMidpoint<dim> q_4;
  deallog << "QMidpoint<" << dim << ">: " << q_4.is_tensor_product()
          << std::endl;

  QSimpson<dim> q_5;
  deallog << "QSimpson<" << dim << ">: " << q_5.is_tensor_product()
          << std::endl;

  QTrapez<dim> q_6;
  deallog << "QTrapez<" << dim << ">: " << q_6.is_tensor_product() << std::endl;

  QMilne<dim> q_7;
  deallog << "QMilne<" << dim << ">: " << q_7.is_tensor_product() << std::endl;

  QWeddle<dim> q_8;
  deallog << "QWeddle<" << dim << ">: " << q_8.is_tensor_product() << std::endl;

  QGaussChebyshev<dim> q_9(3);
  deallog << "QGaussChebyshev<" << dim << ">: " << q_9.is_tensor_product()
          << std::endl;

  QGaussRadauChebyshev<dim> q_10(1);
  deallog << "QGaussRadauChebyshev<" << dim << ">: " << q_10.is_tensor_product()
          << std::endl;

  QGaussLobattoChebyshev<dim> q_11(2);
  deallog << "QGaussLobattoChebyshev<" << dim
          << ">: " << q_11.is_tensor_product() << std::endl;

  QSorted<dim> q_16(QGauss<dim>(3));
  deallog << "QSorted<" << dim << ">: " << q_16.is_tensor_product()
          << std::endl;

  QTelles<dim> q_17(1, Point<dim>());
  deallog << "QTelles<" << dim << ">: " << q_17.is_tensor_product()
          << std::endl;
}

int
main()
{
  initlog();
  deallog << std::boolalpha;

  QGauss<1> q(1);

  print_is_tensor_product<1>();
  QAnisotropic<1> q_1(q);
  deallog << "QAnisotropic<1>: " << q_1.is_tensor_product() << std::endl;
  QGaussLog<1> q_13(1);
  deallog << "QGaussLog<1>: " << q_13.is_tensor_product() << std::endl;
  QGaussLogR<1> q_14(1);
  deallog << "QGaussLogR<1>: " << q_14.is_tensor_product() << std::endl;

  print_is_tensor_product<2>();
  QAnisotropic<2> q_2(q, q);
  deallog << "QAnisotropic<2>: " << q_2.is_tensor_product() << std::endl;
  QGaussOneOverR<2> q_15(2, Point<2>());
  deallog << "QGaussOneOverR<2>: " << q_15.is_tensor_product() << std::endl;

  print_is_tensor_product<3>();
  QAnisotropic<3> q_3(q, q, q);
  deallog << "QAnisotropic<3>: " << q_3.is_tensor_product() << std::endl;
}
