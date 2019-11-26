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

#ifndef dealii_cell_accessor_wrapper_h
#define dealii_cell_accessor_wrapper_h

#include <deal.II/base/config.h>

#include <deal.II/grid/tria_accessor.h>

#include <boost/python.hpp>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class PointWrapper;
  class TriangulationWrapper;

  class CellAccessorWrapper
  {
  public:
    /**
     * Copy constructor.
     */
    CellAccessorWrapper(const CellAccessorWrapper &other);

    /**
     * Constructor. Takes a TriangulationWrapper, the level, and the index
     * associated to the cell.
     */
    CellAccessorWrapper(TriangulationWrapper &triangulation_wrapper,
                        const int             level,
                        const int             index);

    /**
     * Constructor for an empty object.
     */
    CellAccessorWrapper();

    /**
     * Destructor.
     */
    ~CellAccessorWrapper();

    /**
     * Set the refine flag. The possibilities are:
     *  - in 2D: isotropic, no_refinement, cut_x, cut_y, and cut_xy.
     *  - in 3D: isotropic, no_refinement, cut_x, cut_y, cut_z, cut_xy, cut_xz,
     *  cut_yz, and cut_xyz.
     */
    void
    set_refine_flag(const std::string &refinement_case);


    /**
     * Get the refine flag.
     */
    std::string
    get_refine_flag() const;

    /**
     * Set the coarsen flag to true or false.
     */
    void
    set_coarsen_flag(const bool coarsen_flag);

    /**
     * Get the coarsen flag.
     */
    bool
    get_coarsen_flag() const;

    /**
     * Get the barycenter of the cell.
     */
    PointWrapper
    get_barycenter() const;

    /**
     * Set the material id.
     */
    void
    set_material_id(const int material_id);

    /**
     * Get the material id.
     */
    int
    get_material_id() const;

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
     * Return the ith neighbor of a cell. If the neighbor does not exist,
     * i.e., if the ith face of the current object is at the boundary,
     * then an exception is thrown.
     */
    CellAccessorWrapper
    neighbor(const int i) const;

    /**
     * Return whether the cell is at the boundary.
     */
    bool
    at_boundary() const;

    /**
     * This is a slight variation to the at_boundary function:
     * for 2 dimensions it is equivalent, for three
     * dimensions it returns whether at least one of the 12
     * lines of the hexahedron is at a boundary.
     */
    bool
    has_boundary_lines() const;

    /**
     * Get faces of the underlying cell.
     */
    boost::python::list
    faces() const;

    /**
     * Exception.
     */
    DeclException2(ExcVertexDoesNotExist,
                   int,
                   int,
                   << "Requested vertex number " << arg1
                   << " does not exist. The largest vertex number "
                   << "acceptable is " << arg2 - 1);

    DeclException2(ExcNeighborDoesNotExist,
                   int,
                   int,
                   << "Requested neighbor number " << arg1
                   << " does not exist. The largest neighbor number "
                   << "acceptable is " << arg2 - 1);

  private:
    template <int dim_, int spacedim_>
    const CellAccessorWrapper
    construct_neighbor_wrapper(const int i) const;

    /**
     * Dimension of the underlying CellAccessor object.
     */
    int dim;

    /**
     * Space dimension of the underlying CellAccessor object.
     */
    int spacedim;

    /**
     * Pointer that can be casted to the underlying CellAccessor object.
     */
    void *cell_accessor;

    friend class MappingQGenericWrapper;
  };


  template <int dim_, int spacedim_>
  const CellAccessorWrapper
  CellAccessorWrapper::construct_neighbor_wrapper(const int i) const
  {
    Assert(dim_ == dim && spacedim_ == spacedim, ExcInternalError());
    auto *cell = static_cast<CellAccessor<dim_, spacedim_> *>(cell_accessor);

    auto neighbor = cell->neighbor(i);

    CellAccessorWrapper neighbor_wrapper;
    neighbor_wrapper.dim      = dim_;
    neighbor_wrapper.spacedim = spacedim_;
    neighbor_wrapper.cell_accessor =
      new CellAccessor<dim_, spacedim_>(&neighbor->get_triangulation(),
                                        neighbor->level(),
                                        neighbor->index());
    return neighbor_wrapper;
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
