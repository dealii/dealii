// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_tria_accessor_wrapper_h
#define dealii_tria_accessor_wrapper_h

#include <deal.II/base/config.h>

#include <deal.II/grid/tria_accessor.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class PointWrapper;
  class TriangulationWrapper;

  class TriaAccessorWrapper
  {
  public:
    /**
     * Copy constructor.
     */
    TriaAccessorWrapper(const TriaAccessorWrapper &other);

    /**
     * Constructor.
     */
    TriaAccessorWrapper(void *    tria_accessor,
                        const int structdim,
                        const int dim,
                        const int spacedim);

    /**
     * Destructor.
     */
    ~TriaAccessorWrapper();

    /**
     * Get the barycenter of the cell.
     */
    PointWrapper
    get_barycenter() const;

    /**
     * Set the ith vertex of the cell to @p point_wrapper.
     */
    void
    set_vertex(const int i, PointWrapper &point_wrapper);

    /**
     * Return the ith vertex of the cell.
     */
    PointWrapper
    get_vertex(const int i) const;

    /**
     * Set the manifold id.
     */
    void
    set_manifold_id(const int manifold_id);

    /**
     * Get the manifold id.
     */
    int
    get_manifold_id() const;

    /**
     * Set the boundary id.
     */
    void
    set_boundary_id(const int boundary_id);

    /**
     * Get the boundary id.
     */
    int
    get_boundary_id() const;

    /**
     * Set the boundary id for all objects.
     */
    void
    set_all_boundary_ids(const int boundary_id);

    /**
     * Return whether the cell is at the boundary.
     */
    bool
    at_boundary() const;

    /**
     * Compute the dim-dimensional measure of the object.
     * For a dim-dimensional cell in dim-dimensional space,
     * this equals its volume. On the other hand, for a 2d
     * cell in 3d space, or if the current object pointed to
     * is a 2d face of a 3d cell in 3d space, then the function
     * computes the area the object occupies. For a
     * one-dimensional object, return its length.
     */
    double
    measure() const;

    /**
     * Exception.
     */
    DeclException2(ExcVertexDoesNotExist,
                   int,
                   int,
                   << "Requested vertex number " << arg1
                   << " does not exist. The largest vertex number "
                   << "acceptable is " << arg2 - 1);

  private:
    /**
     * Dimension of the underlying TriaAccessor object
     */
    int structdim;

    /**
     * Dimension of the triangulation which TriaAccessor object
     * is linked to.
     */
    int dim;

    /**
     * Space dimension of the underlying TriaAccessor object.
     */
    int spacedim;

    /**
     * Pointer that can be casted to the underlying TriaAccessor object.
     */
    void *tria_accessor;
  };
} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
