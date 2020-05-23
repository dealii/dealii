// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#include <point_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  namespace internal
  {
    template <int dim>
    boost::python::list
    to_list(const void *point)
    {
      const Point<dim> &p = *static_cast<const Point<dim> *>(point);

      boost::python::list p_list;
      for (int d = 0; d < dim; ++d)
        p_list.append(p[d]);

      return p_list;
    }

    template <int dim>
    double
    distance(const Point<dim> &p1, const Point<dim> &p2)
    {
      return p1.distance(p2);
    }



    template <int dim>
    bool
    not_equal(const Point<dim> &p1, const Point<dim> &p2)
    {
      return (p1 != p2);
    }



    template <int dim>
    bool
    equal(const Point<dim> &p1, const Point<dim> &p2)
    {
      return (p1 == p2);
    }



    template <int dim>
    double
    dot_product(const Point<dim> &p1, const Point<dim> &p2)
    {
      return p1 * p2;
    }



    template <int dim>
    Point<dim>
    add_points(const Point<dim> &p1, const Point<dim> &p2)
    {
      return p1 + p2;
    }



    template <int dim>
    Point<dim>
    subtract_point(const Point<dim> &p1, const Point<dim> &p2)
    {
      if (dim == 2)
        return Point<dim>(p1[0] - p2[0], p1[1] - p2[1]);
      else
        return Point<dim>(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
    }



    template <int dim>
    Point<dim>
    opposite_point(const Point<dim> &p)
    {
      return -p;
    }



    template <int dim>
    Point<dim>
    divide_point(const Point<dim> &p, double factor)
    {
      return p / factor;
    }



    template <int dim>
    Point<dim>
    multiply_point(const Point<dim> &p, double factor)
    {
      return p * factor;
    }



    template <int dim>
    void
    add_and_set(Point<dim> &p1, const Point<dim> &p2)
    {
      p1 += p2;
    }



    template <int dim>
    void
    subtract_and_set(Point<dim> &p1, const Point<dim> &p2)
    {
      p1 -= p2;
    }



    template <int dim>
    void
    multiply_and_set(Point<dim> &p, double factor)
    {
      p *= factor;
    }



    template <int dim>
    void
    divide_and_set(Point<dim> &p, double factor)
    {
      p /= factor;
    }
  } // namespace internal

  PointWrapper::PointWrapper()
    : dim(-1)
    , point(nullptr)
  {}


  PointWrapper::PointWrapper(boost::python::list coord)
  {
    dim = boost::python::len(coord);
    if (dim == 2)
      point = new Point<2>(boost::python::extract<double>(coord[0]),
                           boost::python::extract<double>(coord[1]));
    else if (dim == 3)
      point = new Point<3>(boost::python::extract<double>(coord[0]),
                           boost::python::extract<double>(coord[1]),
                           boost::python::extract<double>(coord[2]));
    else
      AssertThrow(
        false,
        ExcMessage(
          "The list of coordinates must contain two or three elements."));
  }



  template <int d>
  PointWrapper::PointWrapper(const Point<d> &p)
  {
    dim = d;
    if (dim == 2)
      point = new Point<2>(p[0], p[1]);

    else
      point = new Point<3>(p[0], p[1], p[2]);
  }



  PointWrapper::PointWrapper(const PointWrapper &other)
  {
    copy(other);
  }



  PointWrapper::~PointWrapper()
  {
    clear();
    dim = -1;
  }



  boost::python::list
  PointWrapper::to_list() const
  {
    if (dim == 2)
      return internal::to_list<2>(point);
    else
      return internal::to_list<3>(point);
  }



  double
  PointWrapper::distance(const PointWrapper &p) const
  {
    AssertThrow(p.get_dim() == dim,
                ExcMessage("The points do not have the same dimension."));

    if (dim == 2)
      return internal::distance(*static_cast<Point<2> *>(point),
                                *static_cast<const Point<2> *>(p.get_point()));
    else
      return internal::distance(*static_cast<Point<3> *>(point),
                                *static_cast<const Point<3> *>(p.get_point()));
  }



  double
  PointWrapper::norm() const
  {
    if (dim == 2)
      return static_cast<Point<2> *>(point)->norm();
    else
      return static_cast<Point<3> *>(point)->norm();
  }



  double
  PointWrapper::norm_square() const
  {
    if (dim == 2)
      return static_cast<Point<2> *>(point)->norm_square();
    else
      return static_cast<Point<3> *>(point)->norm_square();
  }



  PointWrapper &
  PointWrapper::operator=(const PointWrapper &other)
  {
    clear();
    copy(other);

    return *this;
  }



  bool
  PointWrapper::operator!=(const PointWrapper &p) const
  {
    AssertThrow(p.get_dim() == dim,
                ExcMessage("The points do not have the same dimension."));

    if (dim == 2)
      return internal::not_equal(*static_cast<const Point<2> *>(point),
                                 *static_cast<const Point<2> *>(p.get_point()));
    else
      return internal::not_equal(*static_cast<const Point<3> *>(point),
                                 *static_cast<const Point<3> *>(p.get_point()));
  }



  bool
  PointWrapper::operator==(const PointWrapper &p) const
  {
    AssertThrow(p.get_dim() == dim,
                ExcMessage("The points do not have the same dimension."));

    if (dim == 2)
      return internal::equal(*static_cast<const Point<2> *>(point),
                             *static_cast<const Point<2> *>(p.get_point()));
    else
      return internal::equal(*static_cast<const Point<3> *>(point),
                             *static_cast<const Point<3> *>(p.get_point()));
  }



  double PointWrapper::operator*(const PointWrapper &p) const
  {
    AssertThrow(p.get_dim() == dim,
                ExcMessage("The points do not have the same dimension."));

    if (dim == 2)
      return internal::dot_product(*static_cast<const Point<2> *>(point),
                                   *static_cast<const Point<2> *>(
                                     p.get_point()));
    else
      return internal::dot_product(*static_cast<const Point<3> *>(point),
                                   *static_cast<const Point<3> *>(
                                     p.get_point()));
  }



  PointWrapper
  PointWrapper::operator+(const PointWrapper &p) const
  {
    AssertThrow(p.get_dim() == dim,
                ExcMessage("The points do not have the same dimension."));

    if (dim == 2)
      return PointWrapper(
        internal::add_points(*static_cast<const Point<2> *>(point),
                             *static_cast<const Point<2> *>(p.get_point())));
    else
      return PointWrapper(
        internal::add_points(*static_cast<const Point<3> *>(point),
                             *static_cast<const Point<3> *>(p.get_point())));
  }



  PointWrapper
  PointWrapper::operator-(const PointWrapper &p) const
  {
    AssertThrow(p.get_dim() == dim,
                ExcMessage("The points do not have the same dimension."));

    if (dim == 2)
      return PointWrapper(internal::subtract_point(
        *static_cast<const Point<2> *>(point),
        *static_cast<const Point<2> *>(p.get_point())));
    else
      return PointWrapper(internal::subtract_point(
        *static_cast<const Point<3> *>(point),
        *static_cast<const Point<3> *>(p.get_point())));
  }



  PointWrapper
  PointWrapper::operator-() const
  {
    if (dim == 2)
      return PointWrapper(
        internal::opposite_point(*static_cast<Point<2> *>(point)));
    else
      return PointWrapper(
        internal::opposite_point(*static_cast<Point<3> *>(point)));
  }



  PointWrapper
  PointWrapper::operator/(const double factor) const
  {
    AssertThrow(factor != 0., ExcMessage("Dividing by zero."));

    if (dim == 2)
      return PointWrapper(
        internal::divide_point(*static_cast<Point<2> *>(point), factor));
    else
      return PointWrapper(
        internal::divide_point(*static_cast<Point<3> *>(point), factor));
  }


  PointWrapper PointWrapper::operator*(const double factor) const
  {
    if (dim == 2)
      return PointWrapper(
        internal::multiply_point(*static_cast<Point<2> *>(point), factor));
    else
      return PointWrapper(
        internal::multiply_point(*static_cast<Point<3> *>(point), factor));
  }



  PointWrapper &
  PointWrapper::operator+=(const PointWrapper &p)
  {
    AssertThrow(p.get_dim() == dim,
                ExcMessage("The points do not have the same dimension."));

    if (dim == 2)
      internal::add_and_set(*static_cast<Point<2> *>(point),
                            *static_cast<const Point<2> *>(p.get_point()));
    else
      internal::add_and_set(*static_cast<Point<3> *>(point),
                            *static_cast<const Point<3> *>(p.get_point()));

    return *this;
  }



  PointWrapper &
  PointWrapper::operator-=(const PointWrapper &p)
  {
    AssertThrow(p.get_dim() == dim,
                ExcMessage("The points do not have the same dimension."));

    if (dim == 2)
      internal::subtract_and_set(*static_cast<Point<2> *>(point),
                                 *static_cast<const Point<2> *>(p.get_point()));
    else
      internal::subtract_and_set(*static_cast<Point<3> *>(point),
                                 *static_cast<const Point<3> *>(p.get_point()));

    return *this;
  }



  PointWrapper &
  PointWrapper::operator*=(const double factor)
  {
    if (dim == 2)
      internal::multiply_and_set(*static_cast<Point<2> *>(point), factor);
    else
      internal::multiply_and_set(*static_cast<Point<3> *>(point), factor);

    return *this;
  }



  PointWrapper &
  PointWrapper::operator/=(const double factor)
  {
    AssertThrow(factor != 0., ExcMessage("Dividing by zero."));

    if (dim == 2)
      internal::divide_and_set(*static_cast<Point<2> *>(point), factor);
    else
      internal::divide_and_set(*static_cast<Point<3> *>(point), factor);

    return *this;
  }



  double
  PointWrapper::get_x() const
  {
    if (dim == 2)
      return (*static_cast<Point<2> *>(point))(0);
    else
      return (*static_cast<Point<3> *>(point))(0);
  }



  void
  PointWrapper::set_x(double x)
  {
    if (dim == 2)
      (*static_cast<Point<2> *>(point))(0) = x;
    else
      (*static_cast<Point<3> *>(point))(0) = x;
  }



  double
  PointWrapper::get_y() const
  {
    if (dim == 2)
      return (*static_cast<Point<2> *>(point))(1);
    else
      return (*static_cast<Point<3> *>(point))(1);
  }



  void
  PointWrapper::set_y(double y)
  {
    if (dim == 2)
      (*static_cast<Point<2> *>(point))(1) = y;
    else
      (*static_cast<Point<3> *>(point))(1) = y;
  }



  double
  PointWrapper::get_z() const
  {
    if (dim == 3)
      return (*static_cast<Point<3> *>(point))(2);
    else
      AssertThrow(
        false,
        ExcMessage(
          "The z coordinate is only available for three-dimensional points"));
    // Silence a warning
    return 0.;
  }



  void
  PointWrapper::set_z(double z)
  {
    if (dim == 3)
      (*static_cast<Point<3> *>(point))(2) = z;
    else
      AssertThrow(
        false,
        ExcMessage(
          "The z coordinate is only available for three-dimensional points"));
  }



  void
  PointWrapper::clear()
  {
    if (point != nullptr)
      {
        if (dim == 2)
          {
            // We cannot call delete on a void pointer so cast the void pointer
            // back first.
            Point<2> *tmp = static_cast<Point<2> *>(point);
            delete tmp;
          }
        else
          {
            Point<3> *tmp = static_cast<Point<3> *>(point);
            delete tmp;
          }
        point = nullptr;
      }
  }



  void
  PointWrapper::copy(const PointWrapper &other)
  {
    dim = other.dim;

    AssertThrow(other.point != nullptr,
                ExcMessage("Underlying point does not exist."));

    if (dim == 2)
      {
        Point<2> *other_point = static_cast<Point<2> *>(other.point);
        point = new Point<2>((*other_point)[0], (*other_point)[1]);
      }
    else if (dim == 3)
      {
        Point<3> *other_point = static_cast<Point<3> *>(other.point);
        point =
          new Point<3>((*other_point)[0], (*other_point)[1], (*other_point)[2]);
      }
    else
      AssertThrow(false,
                  ExcMessage("The dimension of the point should be 2 or 3."));
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE
