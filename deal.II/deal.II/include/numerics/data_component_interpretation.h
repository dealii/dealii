//---------------------------------------------------------------------------
//    $Id: data_out.h 15313 2007-10-14 01:28:54Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__data_component_interpretation_h
#define __deal2__data_component_interpretation_h



#include <base/config.h>

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
				    * See the @ref step_22 "step-22" tutorial
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
