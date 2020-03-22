/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file ObjectiveFunctionTemplate.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ObjectiveFunctionTemplate.hpp"
#include "QualityMetric.hpp"
#include "MsqError.hpp"
#include "MsqHessian.hpp"
#include "PatchData.hpp"
#include "ElementPatches.hpp"
#include "VertexPatches.hpp"
#include "PatchIterator.hpp"
#include <memory>

namespace MESQUITE_NS {

ObjectiveFunctionTemplate::~ObjectiveFunctionTemplate() {}

void ObjectiveFunctionTemplate::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                                  const Settings* settings,
                                                  MsqError& err )
{
  qualityMetric->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
}

bool ObjectiveFunctionTemplate::initialize_block_coordinate_descent( MeshDomainAssoc* mesh_and_domain,
                                                      const Settings* settings,
                                                      PatchSet* ,
                                                      MsqError& err )
{
  std::auto_ptr<PatchSet> patch_set;
  switch (get_quality_metric()->get_metric_type())
  {
    case QualityMetric::VERTEX_BASED:  
      patch_set = std::auto_ptr<PatchSet>(new VertexPatches( 1, false ));
      break;
    case QualityMetric::ELEMENT_BASED: 
      patch_set = std::auto_ptr<PatchSet>(new ElementPatches);
      break;
    default: 
      MSQ_SETERR(err)("Cannot initialize for BCD for unknown metric type", 
                      MsqError::INVALID_STATE);
      return false;
  }

  Mesh* mesh = mesh_and_domain->get_mesh();
  MeshDomain* domain = mesh_and_domain->get_domain();

  clear();
  patch_set->set_mesh( mesh );
  PatchIterator patches( patch_set.get() );
  
  PatchData pd;
  pd.set_mesh( mesh );
  pd.set_domain( domain );
  if (settings)
    pd.attach_settings( settings );
  
  bool result = true;
  while (patches.get_next_patch( pd, err ) && !MSQ_CHKERR(err))
  {
    double value;
    bool b = evaluate( ObjectiveFunction::ACCUMULATE, pd, value, false, err ); 
    MSQ_ERRZERO(err);
    result = result && b;
  }
  return result;
}

int ObjectiveFunctionTemplate::min_patch_layers() const
{
  if (!get_quality_metric())
    return 0;
  else if (get_quality_metric()->get_metric_type() == QualityMetric::VERTEX_BASED)
    return 2;
  else
    return 1;
}


} // namespace Mesquite
