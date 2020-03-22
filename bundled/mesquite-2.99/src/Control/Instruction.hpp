/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

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

#ifndef MSQ_INSTRUCTION_HPP
#define MSQ_INSTRUCTION_HPP

#include "Mesquite.hpp"

#include <string>

namespace MESQUITE_NS {

class Mesh;
class ParallelMesh;
class MeshDomain;
class MeshDomainAssoc;
class MsqError;
class Settings;

//!\brief Base class for all objects inserted into InstructionQueue
class MESQUITE_EXPORT Instruction
{
  public:
  
    virtual ~Instruction();
    
      //! Called for all instructions in queue before loop_over_mesh
      //! is called for any insetruction in queue.  Default behavior
      //! is to do nothing.
    virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError& err ) = 0;
  
      //! Virtual fuction implementing primary functionaliy of 
      //! instruction instance.
    virtual double loop_over_mesh( MeshDomainAssoc* mesh_and_domain, 
                                   const Settings* settings,
                                   MsqError& err ) = 0;

      //! Virtual fuction implementing primary functionaliy of 
      //! instruction instance for parallel mesh.
    virtual double loop_over_mesh( ParallelMesh* mesh, 
                                   MeshDomain* domain, 
                                   const Settings* settings,
                                   MsqError& err );

      //! Get string name for use in diagnostic and status output
    virtual std::string get_name() const = 0;
    
      //! Using data from query methods in MeshInterface or calculating
      //! as necessary depending on values in Settings, initialize
      //! MSQ_HARD_FIXED and MSQ_SLAVED flags on vertices, and clear
      //! MSQ_CULLED flag on all vertices.  
      //! \NOTE SLAVE_ALL setting is handled in PatchData rather than here 
      //! for efficiency.
    static
    void initialize_vertex_byte( MeshDomainAssoc* mesh_and_domain,
                                 const Settings* settings,
                                 MsqError& err );
};

} // namespace Mesquite

#endif
