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


/** \file IQInterface.hpp
 *  \brief Interface for InstructionQueue and similar classes
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_IQINTERFACE_HPP
#define MSQ_IQINTERFACE_HPP

#include "Mesquite.hpp"
#include "Settings.hpp"
#include "MeshInterface.hpp"

namespace MESQUITE_NS {

class MsqError;
class ParallelMesh;

class MESQUITE_EXPORT IQInterface : public Settings
{
  public:
  
    virtual ~IQInterface();
  
    inline void 
    run_instructions( MeshDomainAssoc* mesh_and_domain, MsqError &err)
      { this->run_common( mesh_and_domain, 0, this, err ); }
    
    inline void 
    run_instructions( Mesh* mesh, MsqError& err )
      { 
        MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(mesh, 0);
        this->run_common( &mesh_and_domain, 0, this, err ); 
      }
    
    inline void 
    run_instructions( ParallelMesh* mesh, MeshDomain* domain, MsqError &err)
      {
        MeshDomainAssoc mesh_and_domain = MeshDomainAssoc((Mesh*)mesh, domain);
        this->run_common( &mesh_and_domain, mesh, this, err ); 
      }
    
    inline void 
    run_instructions( ParallelMesh* mesh, MsqError& err )
      {
        MeshDomainAssoc mesh_and_domain = MeshDomainAssoc((Mesh*)mesh, 0);
        this->run_common( &mesh_and_domain, mesh, this, err ); 
      }

  protected:
  
    virtual void run_common( MeshDomainAssoc* mesh_and_domain,
                             ParallelMesh* pmesh,
                             Settings* settings,
                             MsqError& err ) = 0;
};

} // namespace MESQUITE_NS

#endif
