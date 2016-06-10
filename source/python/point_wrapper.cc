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

#include <deal.II/python/point_wrapper.h>
#include <deal.II/base/point.h>

namespace PyDealII
{
  PointWrapper::PointWrapper()
    :
    dim(-1),
    point(nullptr)
  {}


  PointWrapper::PointWrapper(boost::python::list coord)
  {
    dim = boost::python::len(coord);
    if (dim == 2)
      point = new dealii::Point<2> (boost::python::extract<double>(coord[0]),
                                    boost::python::extract<double>(coord[1]));
    else if (dim == 3)
      point = new dealii::Point<3> (boost::python::extract<double>(coord[0]),
                                    boost::python::extract<double>(coord[1]),
                                    boost::python::extract<double>(coord[2]));
    else
      AssertThrow(false,
                  dealii::ExcMessage("The list of coordinates must contain two or three elements."));
  }



  PointWrapper::~PointWrapper()
  {
    if (point != nullptr)
      {
        if (dim == 2)
          {
            // We cannot call delete on a void pointer so cast the void pointer back
            // first.
            dealii::Point<2> *tmp = static_cast<dealii::Point<2>*>(point);
            delete tmp;
          }
        else
          {
            dealii::Point<3> *tmp = static_cast<dealii::Point<3>*>(point);
            delete tmp;
          }
        point = nullptr;
      }
    dim = -1;
  }



  double PointWrapper::get_x()
  {
    if (dim == 2)
      return (*static_cast<dealii::Point<2>*>(point))(0);
    else
      return (*static_cast<dealii::Point<3>*>(point))(0);
  }



  void PointWrapper::set_x(double x)
  {
    if (dim == 2)
      (*static_cast<dealii::Point<2>*>(point))(0) = x;
    else
      (*static_cast<dealii::Point<3>*>(point))(0) = x;
  }



  double PointWrapper::get_y()
  {
    if (dim == 2)
      return (*static_cast<dealii::Point<2>*>(point))(1);
    else
      return (*static_cast<dealii::Point<3>*>(point))(1);
  }



  void PointWrapper::set_y(double y)
  {
    if (dim == 2)
      (*static_cast<dealii::Point<2>*>(point))(1) = y;
    else
      (*static_cast<dealii::Point<3>*>(point))(1) = y;
  }



  double PointWrapper::get_z()
  {
    if (dim == 3)
      return (*static_cast<dealii::Point<3>*>(point))(2);
    else
      AssertThrow(false,
                  dealii::ExcMessage("The z coordinate is only available for three-dimensional points"));
  }



  void PointWrapper::set_z(double z)
  {
    if (dim == 3)
      (*static_cast<dealii::Point<3>*>(point))(2) = z;
    else
      AssertThrow(false,
                  dealii::ExcMessage("The z coordinate is only available for three-dimensional points"));
  }
}
