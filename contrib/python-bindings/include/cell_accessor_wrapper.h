// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
  class ReferenceCellWrapper;

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

    /*! @copydoc CellAccessor::set_coarsen_flag
     */
    void
    set_coarsen_flag(const bool coarsen_flag);

    /*! @copydoc CellAccessor::coarsen_flag_set
     */
    bool
    get_coarsen_flag() const;

    /*! @copydoc TriaAccessor::barycenter
     */
    PointWrapper
    get_barycenter() const;

    /*! @copydoc TriaAccessor::center
     */
    PointWrapper
    get_center(const bool respect_manifold             = false,
               const bool interpolate_from_surrounding = false) const;

    /*! @copydoc CellAccessor::set_material_id
     */
    void
    set_material_id(const int material_id);

    /*! @copydoc CellAccessor::material_id
     */
    int
    get_material_id() const;

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

    /*! @copydoc TriaAccessor::set_all_manifold_ids
     */
    void
    set_all_manifold_ids(const int manifold_id);

    /*! @copydoc CellAccessor::neighbor
     */
    CellAccessorWrapper
    neighbor(const int i) const;

    /*! @copydoc CellAccessor::at_boundary
     */
    bool
    at_boundary() const;

    /*! @copydoc CellAccessor::has_boundary_lines
     */
    bool
    has_boundary_lines() const;

    /**
     * Get faces of the underlying cell.
     */
    boost::python::list
    faces() const;

    /*! @copydoc TriaAccessor::measure
     */
    double
    measure() const;

    /*! @copydoc CellAccessor::active
     */
    bool
    active() const;

    /*! @copydoc TriaAccessor::level
     */
    int
    level() const;

    /*! @copydoc TriaAccessor::index
     */
    int
    index() const;

    /*! @copydoc CellAccessor::neighbor_is_coarser
     */
    bool
    neighbor_is_coarser(const unsigned int neighbor) const;

    /*! @copydoc CellAccessor::neighbor_of_neighbor
     */
    unsigned int
    neighbor_of_neighbor(const unsigned int neighbor) const;

    /*! @copydoc TriaAccessor::vertex_index
     */
    unsigned int
    vertex_index(const unsigned int i) const;

    /*! @copydoc TriaAccessor::reference_cell
     */
    ReferenceCellWrapper
    reference_cell() const;

    /*! @copydoc TriaAccessor::n_vertices
     */
    unsigned int
    n_vertices() const;

    /*! @copydoc TriaAccessor::n_lines
     */
    unsigned int
    n_lines() const;

    /*! @copydoc TriaAccessor::n_faces
     */
    unsigned int
    n_faces() const;

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

    friend class MappingQWrapper;
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
