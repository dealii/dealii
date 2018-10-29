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
     * Destructor.
     */
    ~CellAccessorWrapper();

    /**
     * Set the refine flag. The possibilities are:
     *  - in 2D: isotropic, no_refinement, cut_x, cut_y, and cut_xy.
     *  - in 3D: isotropic, no_refinement, cut_x, cut_y, cut_z, cut_xy, cut_xz,
     *  cut_yz, and cut_xyz.
     */
    void set_refine_flag(const std::string &refinement_case);


    /**
     * Get the refine flag.
     */
    std::string get_refine_flag() const;

    /**
     * Set the coarsen flag to true or false.
     */
    void set_coarsen_flag(const bool coarsen_flag);

    /**
     * Get the coarsen flag.
     */
    bool get_coarsen_flag() const;

    /**
     * Get the barycenter of the cell.
     */
    PointWrapper get_barycenter() const;

    /**
     * Set the material id.
     */
    void set_material_id(const int material_id);

    /**
     * Get the material id.
     */
    int get_material_id() const;

    /**
     * Set the ith vertex of the cell to @p point_wrapper.
     */
    void set_vertex(const int i, PointWrapper &point_wrapper);

    /**
     * Return the ith vertex of the cell.
     */
    PointWrapper get_vertex(const int i) const;

    /**
     * Set the manifold id.
     */
    void set_manifold_id(const int manifold_id);

    /**
     * Get the manifold id.
     */
    int get_manifold_id() const;

    /**
     * Exception.
     */
    DeclException2(ExcVertexDoesNotExist,
                   int, int,
                   << "Requested vertex number " << arg1
                   << " does not exist. The largest vertex number "
                   << "acceptable is "<< arg2-1);

  private:
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

  };
}

DEAL_II_NAMESPACE_CLOSE

#endif
