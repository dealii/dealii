// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_q1_h
#define dealii_mapping_q1_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup mapping
 * @{
 */

/**
 * Implementation of a $d$-linear mapping from the reference cell to a general
 * quadrilateral/hexahedron.
 *
 * The mapping implemented by this class maps the reference (unit) cell to a
 * general grid cell with straight lines in $d$ dimensions. (Note, however,
 * that in 3d the <i>faces</i> of a general, trilinearly mapped cell may be
 * curved, even if the edges are not). This is the standard mapping used for
 * polyhedral domains. It is also the mapping used throughout deal.II for many
 * functions that come in two variants, one that allows to pass a mapping
 * argument explicitly and one that simply falls back to the MappingQ1 class
 * declared here. (Or, in fact, to an object of kind MappingQ(1), which
 * implements exactly the functionality of this class.)
 *
 * The shape functions for this mapping are the same as for the finite element
 * FE_Q of polynomial degree 1. Therefore, coupling these two yields an
 * isoparametric element.
 *
 * @note This class is, in all reality, nothing more than a different name for
 * calling MappingQ with a polynomial degree of one as argument.
 */
template <int dim, int spacedim = dim>
class MappingQ1 : public MappingQ<dim, spacedim>
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



/**
 * Many places in the library by default use (bi-,tri-)linear mappings unless
 * users explicitly provide a different mapping to use. In these cases, the
 * called function has to create a $Q_1$ mapping object, i.e., an object of
 * kind MappingQ(1). This is costly. It would also be costly to create
 * such objects as static objects in the affected functions, because static
 * objects are never destroyed throughout the lifetime of a program, even
 * though they only have to be created once the first time code runs through a
 * particular function.
 *
 * In order to avoid creation of (static or dynamic) $Q_1$ mapping objects in
 * these contexts throughout the library, this class defines a static $Q_1$
 * mapping object. This object can then be used in all of those places where
 * such an object is needed.
 *
 * @note The use of this object should be avoided since it is only applicable
 *   in cases where a mesh consists exclusively of quadrilaterals or hexahedra.
 *   Use
 * `ReferenceCells::get_hypercube<dim>().get_default_linear_mapping()`
 *   instead.
 */
template <int dim, int spacedim = dim>
struct StaticMappingQ1
{
  /**
   * The static $Q_1$ mapping object discussed in the documentation of this
   * class.
   */
  static MappingQ<dim, spacedim> mapping;
};


/** @} */

template <int dim, int spacedim>
MappingQ<dim, spacedim> StaticMappingQ1<dim, spacedim>::mapping =
  MappingQ1<dim, spacedim>{};

DEAL_II_NAMESPACE_CLOSE

#endif
