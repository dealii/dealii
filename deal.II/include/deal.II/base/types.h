//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2009 by the deal.II authors
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
				    */
  typedef unsigned int subdomain_id_t;

				   /**
				    * A special id for an invalid
				    * subdomain id. This value may not
				    * be used as a valid id but is
				    * used, for example, for default
				    * arguments to indicate a
				    * subdomain id that is not to be
				    * used.
				    *
				    * See the @ref GlossSubdomainId
				    * "glossary" for more information.
				    */
  const unsigned int invalid_subdomain_id = static_cast<subdomain_id_t>(-1);

				   /**
				    * The subdomain id assigned to a
				    * cell whose true subdomain id we
				    * don't know, for example because
				    * it resides on a different
				    * processor on a mesh that is kept
				    * distributed on many
				    * processors. Such cells are
				    * called "artificial".
				    *
				    * See the glossary entries on @ref
				    * GlossSubdomainId "subdomain ids"
				    * and @ref GlossArtificialCell
				    * "artificial cells" as well as
				    * the @ref distributed module for
				    * more information.
				    */
  const unsigned int artificial_subdomain_id = static_cast<subdomain_id_t>(-2);
}


DEAL_II_NAMESPACE_CLOSE

#endif
