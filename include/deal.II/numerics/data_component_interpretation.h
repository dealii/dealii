// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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

#ifndef __deal2__data_component_interpretation_h
#define __deal2__data_component_interpretation_h



#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace solely for the declaration of the
 * DataComponentInterpretation::DataComponentInterpretation enum.
 */
namespace DataComponentInterpretation
{
  /**
   * The members of this enum are used to
   * describe the logical interpretation of
   * what the various components of a
   * vector-valued data set mean. For
   * example, if one has a finite element
   * for the Stokes equations in 2d,
   * representing components $(u,v,p)$, one
   * would like to indicate that the first
   * two, $u$ and $v$, represent a logical
   * vector so that later on when we
   * generate graphical output we can hand
   * them off to a visualization program
   * that will automatically know to render
   * them as a vector field, rather than as
   * two separate and independent scalar
   * fields.
   *
   * By passing a set of enums of the
   * current kind to the
   * DataOut_DoFData::add_data_vector
   * functions, this can be achieved.
   *
   * See the step-22 tutorial
   * program for an example on how this
   * information can be used in practice.
   *
   * @author Wolfgang Bangerth, 2007
   */
  enum DataComponentInterpretation
  {
    /**
     * Indicates that a component of a
     * data set corresponds to a scalar
     * field independent of the others.
     */
    component_is_scalar,

    /**
     * Indicates that a component of a
     * data set is part of a
     * vector-valued quantity.
     */
    component_is_part_of_vector
  };
}


DEAL_II_NAMESPACE_CLOSE

#endif
