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

#ifndef dealii_manifold_wrapper_h
#define dealii_manifold_wrapper_h

#include <deal.II/base/config.h>

#include <boost/python.hpp>

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
     * Create SphericalManifold
     */
    void
    create_spherical(const PointWrapper center);

    /**
     * Create PolarManifold
     */
    void
    create_polar(const PointWrapper center);

    /**
     * Create CylindricalManifold
     */
    void
    create_cylindrical(const int axis, const double tolerance);

    /**
     * Return pointer to an underlying manifold object
     */
    void *
    get_manifold() const;


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
