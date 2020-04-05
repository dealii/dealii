// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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


#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <memory>

#include "../tests.h"


int
main()
{
  initlog();

  int                  a0(-1);
  unsigned int         a1(1);
  double               a2(2.0);
  std::complex<double> a3(3.0, 4.0);
  std::string          a4("foo");
  Point<1>             a5(5.0);
  Point<2>             a6(6.0, 7.0);
  Point<3>             a7(8.0, 9.0, 10.0);

  std::vector<int>                  va0(2, -1);
  std::vector<unsigned int>         va1(2, 1);
  std::vector<double>               va2(2, 2.0);
  std::vector<std::complex<double>> va3(2, std::complex<double>(3.0, 4.0));
  std::vector<std::string>          va4(2, std::string("foo"));
  std::vector<Point<1>>             va5(2, Point<1>(5.0));
  std::vector<Point<2>>             va6(2, Point<2>(6.0, 7.0));
  std::vector<Point<3>>             va7(2, Point<3>(8.0, 9.0, 10.0));



  ParameterHandler prm;
  prm.add_parameter("A signed integer", a0);
  prm.add_parameter("An unsigned integer", a1);
  prm.add_parameter("A double", a2);
  prm.add_parameter("A complex", a3);
  prm.add_parameter("A string", a4);
  prm.add_parameter("A Point<1>", a5);
  prm.add_parameter("A Point<2>", a6);
  prm.add_parameter("A Point<3>", a7);
  prm.add_parameter("A list of signed integers", va0);
  prm.add_parameter("A list of unsigned integers", va1);
  prm.add_parameter("A list of doubles", va2);
  prm.add_parameter("A list of complexes", va3);
  prm.add_parameter("A list of strings", va4);
  prm.add_parameter("A list of Point<1>", va5);
  prm.add_parameter("A list of Point<2>", va6);
  prm.add_parameter("A list of Point<3>", va7);

  prm.log_parameters(deallog);

  prm.set("A signed integer", "-2");
  prm.set("An unsigned integer", "3");
  prm.set("A double", "4.0");
  prm.set("A complex", "5.0, 6.0");
  prm.set("A string", "bar");
  prm.set("A Point<1>", "7.0");
  prm.set("A Point<2>", "8.0, 9.0");
  prm.set("A Point<3>", "10.0, 11.0, 12.0");
  prm.set("A list of signed integers", "-3,-4");
  prm.set("A list of unsigned integers", "6,7");
  prm.set("A list of doubles", "9.0,10.0");
  prm.set("A list of complexes", "1.0,2.0; 3.0,4.0");
  prm.set("A list of strings", "foo, bar");
  prm.set("A list of Point<1>", "1.0; 2.0");
  prm.set("A list of Point<2>", "3.0, 4.0;  5.0,6.0");
  prm.set("A list of Point<3>", "7.0, 8.0, 9.0;  10.0, 11.0, 12.0");

  deallog << "After ParameterHandler::set =========================="
          << std::endl
          << std::endl;
  prm.log_parameters(deallog);

  deallog << "Actual variables            =========================="
          << std::endl
          << std::endl;

  deallog << "Values: " << a0 << std::endl;
  deallog << "Values: " << a1 << std::endl;
  deallog << "Values: " << a2 << std::endl;
  deallog << "Values: " << a3 << std::endl;
  deallog << "Values: " << a4 << std::endl;
  deallog << "Values: " << a5 << std::endl;
  deallog << "Values: " << a6 << std::endl;
  deallog << "Values: " << a7 << std::endl;

  deallog << "Vector values: " << va0[0] << ", " << va0[1] << std::endl;
  deallog << "Vector values: " << va1[0] << ", " << va1[1] << std::endl;
  deallog << "Vector values: " << va2[0] << ", " << va2[1] << std::endl;
  deallog << "Vector values: " << va3[0] << ", " << va3[1] << std::endl;
  deallog << "Vector values: " << va4[0] << ", " << va4[1] << std::endl;
  deallog << "Vector values: " << va5[0] << ", " << va5[1] << std::endl;
  deallog << "Vector values: " << va6[0] << ", " << va6[1] << std::endl;
  deallog << "Vector values: " << va7[0] << ", " << va7[1] << std::endl;

  return 0;
}
