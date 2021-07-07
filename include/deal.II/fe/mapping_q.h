// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2021 by the deal.II authors
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

#ifndef dealii_mapping_q_h
#define dealii_mapping_q_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q_generic.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/

/**
 * A class that implements a polynomial mapping $Q_p$ of degree $p$ on all
 * cells. This class is completely equivalent to the MappingQGeneric class.
 */
template <int dim, int spacedim = dim>
class MappingQ : public MappingQGeneric<dim, spacedim>
{
public:
  /**
   * Constructor.  @p polynomial_degree denotes the polynomial degree of the
   * polynomials that are used to map cells boundary.
   */
  MappingQ(const unsigned int polynomial_degree);

  /**
   * The second argument is here for backward compatibility with previous
   * versions of deal.II, but it does not have any effect on the workings of
   * this class.
   */
  DEAL_II_DEPRECATED_EARLY
  MappingQ(const unsigned int polynomial_degree,
           const bool         use_mapping_q_on_all_cells);

  /**
   * Copy constructor.
   */
  MappingQ(const MappingQ<dim, spacedim> &mapping);
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
