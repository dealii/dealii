// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
