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

#ifndef dealii_mapping_q1_h
#define dealii_mapping_q1_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/


/**
 * Implementation of a $d$-linear mapping from the reference cell to a general
 * quadrilateral/hexahedron.
 *
 * The mapping implemented by this class maps the reference (unit) cell to a
 * general grid cell with straight lines in $d$ dimensions. (Note, however,
 * that in 3D the <i>faces</i> of a general, trilinearly mapped cell may be
 * curved, even if the edges are not). This is the standard mapping used for
 * polyhedral domains. It is also the mapping used throughout deal.II for many
 * functions that come in two variants, one that allows to pass a mapping
 * argument explicitly and one that simply falls back to the MappingQ1 class
 * declared here. (Or, in fact, to an object of kind MappingQGeneric(1), which
 * implements exactly the functionality of this class.)
 *
 * The shape functions for this mapping are the same as for the finite element
 * FE_Q of polynomial degree 1. Therefore, coupling these two yields an
 * isoparametric element.
 *
 * @note This class is, in all reality, nothing more than a different name for
 * calling MappingQGeneric with a polynomial degree of one as argument.
 */
template <int dim, int spacedim = dim>
class MappingQ1 : public MappingQGeneric<dim, spacedim>
{
public:
  /**
   * Default constructor.
   */
  MappingQ1();

  // for documentation, see the Mapping base class
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;
};


/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
