// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__mapping_q1_h
#define dealii__mapping_q1_h


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
 *
 * @author Guido Kanschat, 2000, 2001; Ralf Hartmann, 2000, 2001, 2005,
 * Wolfgang Bangerth, 2015
 */
template <int dim, int spacedim=dim>
class MappingQ1 : public MappingQGeneric<dim,spacedim>
{
public:
  /**
   * Default constructor.
   */
  MappingQ1 ();

  // for documentation, see the Mapping base class
  virtual
  MappingQ1<dim,spacedim> *clone () const;
};



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
template <int dim, int spacedim=dim>
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
