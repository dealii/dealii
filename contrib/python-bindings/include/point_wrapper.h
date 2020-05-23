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

#ifndef dealii_point_wrapper_h
#define dealii_point_wrapper_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <boost/python.hpp>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class PointWrapper
  {
  public:
    /**
     * This constructor exists only so that a std::vector of PointWrapper can be
     * created. This constructor should not be used.
     */
    PointWrapper();

    /**
     * Constructor. Use a list of double to initialize the Point.
     */
    PointWrapper(boost::python::list list);

    /**
     * Constructor. Initialize PointWrapper using a Point.
     */
    template <int dim>
    PointWrapper(const Point<dim> &p);

    /**
     * Copy constructor.
     */
    PointWrapper(const PointWrapper &other);

    /**
     * Destructor.
     */
    ~PointWrapper();

    /*! @copydoc Point::distance
     */
    double
    distance(const PointWrapper &p) const;

    /*! @copydoc Tensor::norm
     */
    double
    norm() const;

    /*! @copydoc Tensor::norm_square
     */
    double
    norm_square() const;

    /**
     * Convert point's coordinates to a python list with [x,y] or [x,y,z]
     * for 2-D and 3-D, respectively.
     */
    boost::python::list
    to_list() const;

    /**
     * Assignment operator. The dimension of the point is changed if it is
     * different than the one of @p other.
     */
    PointWrapper &
    operator=(const PointWrapper &other);

    /**
     * Test for inequality of two points.
     */
    bool
    operator!=(const PointWrapper &p) const;

    /**
     * Test for equality of two points.
     */
    bool
    operator==(const PointWrapper &p) const;

    /**
     * Return the scalar product of the vectors representing two points.
     */
    double operator*(const PointWrapper &p) const;

    /**
     * Add an offset to a point.
     */
    PointWrapper
    operator+(const PointWrapper &) const;

    /**
     * Subtract two points.
     */
    PointWrapper
    operator-(const PointWrapper &) const;

    /**
     * The opposite point.
     */
    PointWrapper
    operator-() const;

    /**
     * Divide the coordinates of the point by a factor.
     */
    PointWrapper
    operator/(const double factor) const;

    /**
     * Multiply the coordinates of the point by a factor.
     */
    PointWrapper operator*(const double factor) const;

    /**
     * Add another point.
     */
    PointWrapper &
    operator+=(const PointWrapper &p);

    /**
     * Subtract another point.
     */
    PointWrapper &
    operator-=(const PointWrapper &p);

    /**
     * Scale the coordinates of the point by factor.
     */
    PointWrapper &
    operator*=(const double factor);

    /**
     * Scale the coordinates of the point by 1/factor.
     */
    PointWrapper &
    operator/=(const double factor);

    /**
     * Return the first component of the Point.
     */
    double
    get_x() const;

    /**
     * Set the first component of the Point.
     */
    void
    set_x(double x);

    /**
     * Return the second component of the Point.
     */
    double
    get_y() const;

    /**
     * Set the second component of the Point.
     */
    void
    set_y(double y);

    /**
     * Return the third component of the Point.
     */
    double
    get_z() const;

    /**
     * Set the third component of the Point.
     */
    void
    set_z(double z);

    /**
     * Return the dimension of the underlying Point.
     */
    int
    get_dim() const;

    /**
     * Return a pointer that can be casted to the underlying Point.
     */
    void *
    get_point();

    /**
     * Return a constant pointer that can be casted to the underlying Point.
     */
    const void *
    get_point() const;

  private:
    /**
     * Delete the underlying Point and free the memory.
     */
    void
    clear();


    /**
     * Copy @p other PointWrapper.
     */
    void
    copy(const PointWrapper &other);

    /**
     * Dimension of the Point.
     */
    int dim;

    /**
     * Pointer to the underlying Point object.
     */
    void *point;

    friend class MappingQGenericWrapper;
  };


  //--------------------- Inline functions ----------------------//



  inline int
  PointWrapper::get_dim() const
  {
    return dim;
  }



  inline void *
  PointWrapper::get_point()
  {
    return point;
  }



  inline const void *
  PointWrapper::get_point() const
  {
    return point;
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
