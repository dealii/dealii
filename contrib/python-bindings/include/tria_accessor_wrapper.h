// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

    /*! @copydoc TriaAccessor::barycenter
     */
    PointWrapper
    get_barycenter() const;

    /*! @copydoc TriaAccessor::center
     */
    PointWrapper
    get_center(const bool respect_manifold             = false,
               const bool interpolate_from_surrounding = false) const;

    /**
     * Set the ith vertex of the cell to @p point_wrapper.
     */
    void
    set_vertex(const int i, PointWrapper &point_wrapper);

    /*! @copydoc TriaAccessor::vertex
     */
    PointWrapper
    get_vertex(const int i) const;

    /*! @copydoc TriaAccessor::set_manifold_id
     */
    void
    set_manifold_id(const int manifold_id);

    /*! @copydoc TriaAccessor::manifold_id
     */
    int
    get_manifold_id() const;

    /*! @copydoc TriaAccessor::set_boundary_id
     */
    void
    set_boundary_id(const int boundary_id);

    /*! @copydoc TriaAccessor::boundary_id
     */
    int
    get_boundary_id() const;

    /*! @copydoc TriaAccessor::set_all_boundary_ids
     */
    void
    set_all_boundary_ids(const int boundary_id);

    /*! @copydoc TriaAccessor::at_boundary
     */
    bool
    at_boundary() const;

    /*! @copydoc TriaAccessor::measure
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
