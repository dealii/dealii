// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

#include <boost/python.hpp>
#include <deal.II/base/point.h>

namespace PyDealII
{

  double get_x(const dealii::Point<2> &point)
  {
    return point(0);
  }

  double get_y(const dealii::Point<2> &point)
  {
    return point(1);
  }

  void export_point()
  {
    boost::python::class_<dealii::Point<2>> ("Point", "Constructor. Requires x and y.",
                                             boost::python::init<double, double>())
                                         .add_property("x", get_x, "Return the x component of the point.")
                                         .add_property("y", get_y, "Return the y component of the point.");
  }

}
