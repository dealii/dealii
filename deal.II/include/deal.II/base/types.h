//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2009, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__types_h
#define __deal2__types_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which we declare typedefs for types used in deal.II,
 * as well as special values for these types.
 */
namespace types
{
                                   /**
                                    * The type used to denote
                                    * subdomain_ids of cells.
                                    *
                                    * See the @ref GlossSubdomainId
                                    * "glossary" for more information.
				    *
				    * There is a special value,
				    * numbers::invalid_subdomain_id
				    * that is used to indicate an
				    * invalid value of this type.
                                    */
  typedef unsigned int subdomain_id;

				   /**
				    * @deprecated Old name for the typedef above.
				    */
  typedef subdomain_id subdomain_id_t;

                                   /**
				    * @deprecated Use numbers::invalid_subdomain_id
                                    */
  const unsigned int invalid_subdomain_id = static_cast<subdomain_id>(-1);

                                   /**
                                    * @deprecated Use numbers::artificial_subdomain_id
                                    */
  const unsigned int artificial_subdomain_id = static_cast<subdomain_id>(-2);

                                   /**
                                    * The type used to denote global dof
                                    * indices.
                                    */
  typedef unsigned int global_dof_index;

                                   /**
				    *  @deprecated Use numbers::invalid_dof_index
                                    */
  const global_dof_index invalid_dof_index = static_cast<global_dof_index>(-1);

          /**
           * The type used to denote boundary indicators associated with every
           * piece of the boundary and, in the case of meshes that describe
           * manifolds in higher dimensions, associated with every cell.
	   *
	   * There is a special value, numbers::internal_face_boundary_id
	   * that is used to indicate an invalid value of this type and that
	   * is used as the boundary indicator for faces that are in the interior
	   * of the domain and therefore not part of any addressable boundary
	   * component.
           */
  typedef unsigned char boundary_id;

				   /**
				    * @deprecated Old name for the typedef above.
				    */
  typedef boundary_id boundary_id_t;

          /**
           * The type used to denote material indicators associated with every
           * cell.
	   *
	   * There is a special value, numbers::invalid_material_id
	   * that is used to indicate an invalid value of this type.
           */
  typedef unsigned char material_id;

				   /**
				    * @deprecated Old name for the typedef above.
				    */
  typedef material_id material_id_t;

}


DEAL_II_NAMESPACE_CLOSE

#endif
