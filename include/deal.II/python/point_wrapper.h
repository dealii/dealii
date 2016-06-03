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

#ifndef dealii__point_wrapper_h
#define dealii__point_wrapper_h
#include <boost/python.hpp>

namespace PyDealII
{
  class PointWrapper
  {
  public:
    PointWrapper(boost::python::list list);
    double get_x();
    void set_x(double x);
    double get_y();
    void set_y(double y);
    double get_z();
    void set_z(double z);

  private:
    int dim;
    std::shared_ptr<void> point;
  };
}

#endif
