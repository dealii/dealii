//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__integrators_elasticity_h
#define __deal2__integrators_elasticity_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/quadrature.h>
#include <lac/full_matrix.h>
#include <fe/mapping.h>
#include <fe/fe_values.h>
#include <numerics/mesh_worker_info.h>

DEAL_II_NAMESPACE_OPEN

namespace LocalIntegrators
{
/**
 * @brief Local integrators related to elasticity problems.
 *
 * @ingroup Integrators
 * @author Guido Kanschat
 * @date 2010
 */
  namespace Elasticity
  {
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
