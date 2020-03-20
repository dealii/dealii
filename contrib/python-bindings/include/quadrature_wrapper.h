// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#ifndef dealii_quadrature_wrapper_h
#define dealii_quadrature_wrapper_h

#include <deal.II/base/config.h>

#include <boost/python.hpp>

#include <point_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class QuadratureWrapper
  {
  public:
    /**
     * Copy constructor.
     */
    QuadratureWrapper(const QuadratureWrapper &other);

    /**
     * Constructor.
     */
    QuadratureWrapper(const int dim);

    /**
     * Destructor.
     */
    ~QuadratureWrapper();

    /**
     * Create QGauss.
     */
    void
    create_gauss(const unsigned int n);

    /**
     * Create QGaussLobatto.
     */
    void
    create_gauss_lobatto(const unsigned int n);

    /*! @copydoc Quadrature::get_points
     */
    boost::python::list
    get_points() const;

    /*! @copydoc Quadrature::get_weights
     */
    boost::python::list
    get_weights() const;

    /**
     * Return pointer to an underlying Quadrature object
     */
    void *
    get_quadrature() const;

    /**
     * Return the dimension of the underlying object
     */
    int
    get_dim() const;

  private:
    /**
     * Dimension of the underlying Quadrature object.
     */
    int dim;

    /**
     * Pointer that can be casted to the underlying Quadrature object.
     */
    void *quadrature_ptr;
  };

} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
