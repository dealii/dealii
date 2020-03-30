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

#ifndef dealii_mapping_wrapper_h
#define dealii_mapping_wrapper_h

#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <boost/python.hpp>

#include <cell_accessor_wrapper.h>
#include <point_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class MappingQGenericWrapper
  {
  public:
    /**
     * Copy constructor.
     */
    MappingQGenericWrapper(const MappingQGenericWrapper &other);

    /**
     * Default constructor.
     */
    MappingQGenericWrapper();

    /**
     * Constructor.
     */
    MappingQGenericWrapper(const int dim, const int spacedim, const int degree);

    /**
     * Destructor.
     */
    ~MappingQGenericWrapper();

    /*! @copydoc Mapping::transform_unit_to_real_cell
     */
    PointWrapper
    transform_unit_to_real_cell(CellAccessorWrapper &cell, PointWrapper &point);

    /*! @copydoc Mapping::transform_real_to_unit_cell
     */
    PointWrapper
    transform_real_to_unit_cell(CellAccessorWrapper &cell, PointWrapper &point);

    /*! @copydoc Mapping::project_real_point_to_unit_point_on_face
     */
    PointWrapper
    project_real_point_to_unit_point_on_face(CellAccessorWrapper &cell,
                                             const unsigned int   face_no,
                                             PointWrapper &       point);

    /**
     * Get the underlying mapping.
     */
    void *
    get_mapping() const;

  private:
    /**
     * Dimension of the underlying Mapping object.
     */
    int dim;

    /**
     * Space dimension of the underlying Mapping object.
     */
    int spacedim;

    /**
     * Degree of the underlying mapping.
     */
    int degree;

    /**
     * Pointer that can be casted to the underlying Mapping object.
     */
    void *mapping_ptr;
  };

} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
