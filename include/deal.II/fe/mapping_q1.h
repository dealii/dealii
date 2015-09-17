// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/


/**
 * Implementation of a $d$-linear mapping from the reference cell to a general
 * quadrilateral/hexahedron.
 *
 * The mapping implemented by this class maps the reference (unit) cell
 * to a general grid cell with
 * straight lines in $d$ dimensions. (Note, however, that in 3D the
 * <i>faces</i> of a general, trilinearly mapped cell may be curved, even if the
 * edges are not). This is the standard mapping used for polyhedral domains. It
 * is also the mapping used throughout deal.II for many functions that come in
 * two variants, one that allows to pass a mapping argument explicitly and one
 * that simply falls back to the MappingQ1 class declared here.
 *
 * The shape functions for this mapping are the same as for the finite element FE_Q
 * of polynomial degree 1. Therefore, coupling these two yields an isoparametric element.
 *
 * @author Guido Kanschat, 2000, 2001; Ralf Hartmann, 2000, 2001, 2005, Wolfgang Bangerth, 2015
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

  /**
   * @name Mapping points between reference and real cells
   * @{
   */

  // for documentation, see the Mapping base class
  virtual
  Point<dim>
  transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<spacedim>                            &p) const;

  /**
   * @}
   */
};



/**
 * In order to avoid creation of static MappingQ1 objects at several places in
 * the library, we define a static MappingQ1 object once and for all, for use
 * in places where a $Q_1$ mapping is required but do not want to create a new
 * object of this type everytime we get there.
 */
template <int dim, int spacedim=dim>
struct StaticMappingQ1
{
  static MappingQ1<dim, spacedim> mapping;
};


/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
