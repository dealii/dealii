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
#include <deal.II/python/point_wrapper.h>
#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>

namespace PyDealII
{
  void export_point()
  {
    boost::python::class_<PointWrapper> ("Point", boost::python::init<boost::python::list>())
    .add_property("x", &PointWrapper::get_x, &PointWrapper::set_x, "Get the x component of the point.")
    .add_property("y", &PointWrapper::get_y, &PointWrapper::set_y, "Get the y component of the point.")
    .add_property("z", &PointWrapper::get_z, &PointWrapper::set_z, "Get the z component of the point.");
  }
}
