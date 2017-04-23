// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2015 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>

#include <memory>

using namespace dealii;

int main()
{
  initlog();

  // create some List patterns and try to convert to and from string
  Patterns::List p0(Patterns::Integer() ,0,10,",");
  Patterns::List p1(Patterns::Integer(0),0,10,",");
  Patterns::List p2(Patterns::Double()  ,0,10,",");
  Patterns::List p3(Patterns::Anything(),0,10,",");
  Patterns::List p4(Patterns::Double(), 3,3  ,","); // A Point<3>

  Patterns::List vp0(Patterns::List(Patterns::Integer()     ),0,10,";");
  Patterns::List vp1(Patterns::List(Patterns::Integer(0)    ),0,10,";");
  Patterns::List vp2(Patterns::List(Patterns::Double()      ),0,10,";");
  Patterns::List vp3(Patterns::List(Patterns::Anything()    ),0,10,";");
  Patterns::List vp4(Patterns::List(Patterns::Double(), 3,3 ),0,10,";");


  std::vector<int>          l0 (2,-1);
  std::vector<unsigned int> l1 (3,2);
  std::vector<double>       l2 (4, 3.14);
  std::vector<std::string>  l3 (2, "bar");
  Point<3>                  l4 (3, 2, 1);

  std::vector<std::vector<int>          > vl0 (2,std::vector<int>          (3,1));
  std::vector<std::vector<unsigned int> > vl1 (2,std::vector<unsigned int> (3,2));
  std::vector<std::vector<double>       > vl2 (2,std::vector<double>       (3, 3.14));
  std::vector<std::vector<std::string>  > vl3 (2,std::vector<std::string>  (3, "foo"));
  std::vector<Point<3>                  > vl4 (2,Point<3>                  (4,3,2));

  deallog << "List of int         : " << p0.to_string(l0) << std::endl;
  deallog << "List of unsigned int: " << p1.to_string(l1) << std::endl;
  deallog << "List of double      : " << p2.to_string(l2) << std::endl;
  deallog << "List of string      : " << p3.to_string(l3) << std::endl;
  deallog << "Point<3>            : " << p4.to_string(l4) << std::endl;

  deallog << "List of lists of int         : " << vp0.to_string(vl0) << std::endl;
  deallog << "List of lists of unsigned int: " << vp1.to_string(vl1) << std::endl;
  deallog << "List of lists of double      : " << vp2.to_string(vl2) << std::endl;
  deallog << "List of lists of string      : " << vp3.to_string(vl3) << std::endl;
  deallog << "List of Point<3>             : " << vp4.to_string(vl4) << std::endl;

  deallog << "=============================" << std::endl;

  p0.to_value("1,2,3"  , l0);
  p1.to_value("3,4,5"  , l1);
  p2.to_value("5,6,7"  , l2);
  p3.to_value("8,9,8.5", l3);
  p4.to_value("8,9,8.5", l4);

  vp0.to_value("1,2,3  ; 1,2,3  ", vl0);
  vp1.to_value("3,4,5  ; 3,4,5  ", vl1);
  vp2.to_value("5,6,7  ; 5,6,7  ", vl2);
  vp3.to_value("8,9,8.5; 8,9,8.5", vl3);
  vp4.to_value("8,9,8.5; 8,9,8.5", vl4);

  deallog << "List of int         : " << p0.to_string(l0) << std::endl;
  deallog << "List of unsigned int: " << p1.to_string(l1) << std::endl;
  deallog << "List of double      : " << p2.to_string(l2) << std::endl;
  deallog << "List of string      : " << p3.to_string(l3) << std::endl;
  deallog << "Point<3>            : " << p4.to_string(l4) << std::endl;

  deallog << "List of lists of int         : " << vp0.to_string(vl0) << std::endl;
  deallog << "List of lists of unsigned int: " << vp1.to_string(vl1) << std::endl;
  deallog << "List of lists of double      : " << vp2.to_string(vl2) << std::endl;
  deallog << "List of lists of string      : " << vp3.to_string(vl3) << std::endl;
  deallog << "List of Point<3>             : " << vp4.to_string(vl4) << std::endl;

  return 0;
}
