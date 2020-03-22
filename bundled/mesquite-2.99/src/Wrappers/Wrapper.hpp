/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Wrapper.hpp
 *  \brief Common interface implemented by wrappers.
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_WRAPPER_HPP
#define MSQ_WRAPPER_HPP

#include "Mesquite.hpp"
#include "IQInterface.hpp"

namespace MESQUITE_NS {

class QualityAssessor;

/**\brief Interface for wrappers.  
 *
 * Interaface implemented by wrappers.  In addition to implementing
 * IQInterface, also provide access to QualityAssessor instance so 
 * that caller can modify QA output.
 */
class MESQUITE_EXPORT Wrapper : public IQInterface
{
  public:
    
    Wrapper();  
  
    virtual ~Wrapper();
    
    /** Get the quality assessor associated with this wrapper */
    inline
    QualityAssessor& quality_assessor()
      { return *qualAssessor; }
    
    /** Get the quality assessor associated with this wrapper */
    inline
    const QualityAssessor& quality_asssessor() const
      { return *qualAssessor; }

  protected:
  
    /** Function inherited from IQInterface that we implement here */
    void run_common( MeshDomainAssoc* mesh_and_domain,
                     ParallelMesh* pmesh,
                     Settings* settings,
                     MsqError& err );
  
    /** Function that each wrapper must implement */
    virtual void run_wrapper( MeshDomainAssoc* mesh_and_domain,
                              ParallelMesh* pmesh,
                              Settings* settings,
                              QualityAssessor* quality_assessor,
                              MsqError& err ) = 0;
    

  private:
    
    QualityAssessor* qualAssessor;
};

} // namespace MESQUITE_NS

#endif
