// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_manifold_wrapper_h
#define dealii_manifold_wrapper_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>

#include <deal.II/lac/vector.h>

#include <boost/python.hpp>

#include <function_wrapper.h>
#include <point_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class ManifoldWrapper
  {
  public:
    /**
     * Copy constructor.
     */
    ManifoldWrapper(const ManifoldWrapper &other);

    /**
     * Constructor.
     */
    ManifoldWrapper(const int dim, const int spacedim);

    /**
     * Destructor.
     */
    ~ManifoldWrapper();

    /**
     * Create SphericalManifold.
     */
    void
    create_spherical(const PointWrapper center);

    /**
     * Create PolarManifold.
     */
    void
    create_polar(const PointWrapper center);

    /**
     * Create CylindricalManifold along the fixed axis.
     */
    void
    create_cylindrical(const int axis = 0, const double tolerance = 1e-10);

    /**
     * Create CylindricalManifold with the given orientation.
     */
    void
    create_cylindrical(const boost::python::list &direction,
                       const boost::python::list &axial_point);

    /**
     * Create FunctionManifold with string expressions for the push
     * forward and pull back functions.
     */
    void
    create_function_string(const std::string &push_forward,
                           const std::string &pull_back);

    /**
     * Create FunctionManifold with python the push forward and
     * pull back functions.
     */
    void
    create_function(boost::python::object &push_forward,
                    boost::python::object &pull_back);

    /**
     * Return pointer to an underlying manifold object
     */
    void *
    get_manifold() const;

    /**
     * Return the dimension of the underlying object
     */
    int
    get_dim() const;

    /**
     * Return the space dimension of the underlying object
     */
    int
    get_spacedim() const;

  private:
    /**
     * Dimension of the underlying Manifold object.
     */
    int dim;

    /**
     * Space dimension of the underlying Manifold object.
     */
    int spacedim;

    /**
     * Pointer that can be casted to the underlying Mapping object.
     */
    void *manifold_ptr;
  };

} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
