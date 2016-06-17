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
     * Copy constructor.
     */
    PointWrapper(const PointWrapper &other);

    /**
     * Destructor.
     */
    ~PointWrapper();

    /**
     * Assignment operator. The dimension of the point is changed if it is
     * different than the one of @p other.
     */
    PointWrapper &operator= (const PointWrapper &other);

    /**
     * Return the first component of the Point.
     */
    double get_x() const;

    /**
     * Set the first component of the Point.
     */
    void set_x(double x);

    /**
     * Return the second component of the Point.
     */
    double get_y() const;

    /**
     * Set the second component of the Point.
     */
    void set_y(double y);

    /**
     * Return the third component of the Point.
     */
    double get_z() const;

    /**
     * Set the third component of the Point.
     */
    void set_z(double z);

    /**
     * Return the dimension of the underlying Point.
     */
    int get_dim() const;

    /**
     * Return a pointer that can be casted to the underlying Point.
     */
    void *get_point();

  private:
    /**
     * Delete the underlying Point and free the memory.
     */
    void clear();


    /**
     * Copy @p other PointWrapper.
     */
    void copy(const PointWrapper &other);

    /**
     * Dimension of the Point.
     */
    int dim;

    /**
     * Pointer to the underlying Point object.
     */
    void *point;
  };


  //--------------------- Inline functions ----------------------//



  inline
  int PointWrapper::get_dim() const
  {
    return dim;
  }



  inline
  void *PointWrapper::get_point()
  {
    return point;
  }
}

#endif
