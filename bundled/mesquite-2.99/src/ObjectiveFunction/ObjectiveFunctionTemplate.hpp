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


/** \file ObjectiveFunctionTemplate.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_OBJECTIVE_FUNCTION_TEMPLATE_HPP
#define MSQ_OBJECTIVE_FUNCTION_TEMPLATE_HPP

#include "Mesquite.hpp"
#include "ObjectiveFunction.hpp"

namespace MESQUITE_NS {

/**\brief Base for most concrete objective functions
 *
 * Base class for objective functions which are a function
 * of a quality metric.
 */
class MESQUITE_EXPORT ObjectiveFunctionTemplate : public ObjectiveFunction
{
  public:
  
    ObjectiveFunctionTemplate( QualityMetric* qm = 0 ) : qualityMetric(qm) {}
    
    virtual ~ObjectiveFunctionTemplate();
    
    QualityMetric* get_quality_metric() const { return qualityMetric; }
    
    void set_quality_metric( QualityMetric* metric ) { qualityMetric = metric; }

    virtual bool initialize_block_coordinate_descent( MeshDomainAssoc* mesh_and_domain,
                                                      const Settings* settings,
                                                      PatchSet* user_set,
                                                      MsqError& err );

    virtual int min_patch_layers() const;
     
      //!\brief Called at start of instruction queue processing
      //!
      //! Do any preliminary global initialization, consistency checking,
      //! etc.  Default implementation does nothing.
     virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                    const Settings* settings,
                                    MsqError& err );

  private:
  
    QualityMetric* qualityMetric;
};

} // namespace Mesquite

#endif
