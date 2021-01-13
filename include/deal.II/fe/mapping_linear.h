// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
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

#ifndef dealii_mapping_linear_h
#define dealii_mapping_linear_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q1.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/


/**
 * Many places in the library by default use (bi-,tri-)linear mappings unless
 * users explicitly provide a different mapping to use. In these cases, the
 * called function has to create a $Q_1$ mapping object, i.e., an object of
 * kind MappingQGeneric(1). This is costly. It would also be costly to create
 * such objects as static objects in the affected functions, because static
 * objects are never destroyed throughout the lifetime of a program, even
 * though they only have to be created once the first time code runs through a
 * particular function.
 *
 * In order to avoid creation of (static or dynamic) $Q_1$ mapping objects in
 * these contexts throughout the library, this class defines a static $Q_1$
 * mapping object. This object can then be used in all of those places where
 * such an object is needed.
 */
template <int dim, int spacedim = dim>
struct StaticMappingQ1
{
  /**
   * The static $Q_1$ mapping object discussed in the documentation of this
   * class.
   */
  static MappingQGeneric<dim, spacedim> mapping;
};


/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
