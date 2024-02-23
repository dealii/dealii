// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_wrapper_h
#define dealii_mapping_wrapper_h

#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q.h>

#include <boost/python.hpp>

#include <cell_accessor_wrapper.h>
#include <point_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class MappingQWrapper
  {
  public:
    /**
     * Copy constructor.
     */
    MappingQWrapper(const MappingQWrapper &other);

    /**
     * Default constructor.
     */
    MappingQWrapper();

    /**
     * Constructor.
     */
    MappingQWrapper(const int dim, const int spacedim, const int degree);

    /**
     * Destructor.
     */
    ~MappingQWrapper();

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
                                             PointWrapper        &point);

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
