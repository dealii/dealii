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
#include <deal.II/base/exceptions.h>

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



  PointWrapper::PointWrapper(boost::python::list coord)
  {
    dim = boost::python::len(coord);
    if (dim == 2)
      point.reset(new dealii::Point<2> (boost::python::extract<double>(coord[0]),
                                        boost::python::extract<double>(coord[1])));
    else if (dim == 3)
      point.reset(new dealii::Point<3> (boost::python::extract<double>(coord[0]),
                                        boost::python::extract<double>(coord[1]),
                                        boost::python::extract<double>(coord[2])));
    else
      AssertThrow(false,
                  dealii::ExcMessage("The list of coordinates must contain two or three elements."));
  }



  double PointWrapper::get_x()
  {
    if (dim == 2)
      return (*std::static_pointer_cast<dealii::Point<2>>(point))(0);
    else
      return (*std::static_pointer_cast<dealii::Point<3>>(point))(0);
  }



  void PointWrapper::set_x(double x)
  {
    if (dim == 2)
      (*std::static_pointer_cast<dealii::Point<2>>(point))(0) = x;
    else
      (*std::static_pointer_cast<dealii::Point<3>>(point))(0) = x;
  }



  double PointWrapper::get_y()
  {
    if (dim == 2)
      return (*std::static_pointer_cast<dealii::Point<2>>(point))(1);
    else
      return (*std::static_pointer_cast<dealii::Point<3>>(point))(1);
  }



  void PointWrapper::set_y(double y)
  {
    if (dim == 2)
      (*std::static_pointer_cast<dealii::Point<2>>(point))(1) = y;
    else
      (*std::static_pointer_cast<dealii::Point<3>>(point))(1) = y;
  }



  double PointWrapper::get_z()
  {
    if (dim == 3)
      return (*std::static_pointer_cast<dealii::Point<3>>(point))(2);
    else
      AssertThrow(false,
                  dealii::ExcMessage("The z coordinate is only available for three-dimensional points"));
  }



  void PointWrapper::set_z(double z)
  {
    if (dim == 3)
      (*std::static_pointer_cast<dealii::Point<3>>(point))(2) = z;
    else
      AssertThrow(false,
                  dealii::ExcMessage("The z coordinate is only available for three-dimensional points"));
  }



  void export_point()
  {
    boost::python::class_<PointWrapper> ("Point",boost::python::init<boost::python::list>())
    .add_property("x", &PointWrapper::get_x, &PointWrapper::set_x, "Get the x component of the point.")
    .add_property("y", &PointWrapper::get_y, &PointWrapper::set_y, "Get the y component of the point.")
    .add_property("z", &PointWrapper::get_z, &PointWrapper::set_z, "Get the z component of the point.");
  }

}
