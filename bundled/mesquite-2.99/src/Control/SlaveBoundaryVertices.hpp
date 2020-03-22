/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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


/** \file SlaveBoundaryVertices.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_SLAVE_BOUNDARY_VERTICES_HPP
#define MSQ_SLAVE_BOUNDARY_VERTICES_HPP

#include "Mesquite.hpp"
#include "VertexSlaver.hpp"

namespace MESQUITE_NS {

/**\brief Utility to set slaved vs. non-slaved vertices.
 *
 * This class can be inserted in the instruction queue before any
 * optimization to determine which higher-order nodes are slaved
 * as a function of their distance from the bounary of the mesh.
 * 
 * This implementation assumes that the boundary of the mesh can
 * be identified either by dimension of the associated geometric
 * domain or by fixed vertices.
 */
class MESQUITE_EXPORT SlaveBoundaryVertices : public VertexSlaver
{
  public:
    /**\brief Define defintion of boundary as either the distance
     *        from a fixed vertex or the distance from the domain
     *        of the boundary vertices.
     *
     * If the second parameter is not specified, the boundary will
     * be assumed to be indicated by fixed vertices.  If it is specified,
     * the boundary will be assumed to be those vertices constrained to
     * a domain less than or equal to the specified dimension (I.e. 
     * vertices constrained such that the number of degrees of freedom
     * in their motion is less than or equal to the specified number.)
     *
     *\param depth The number of elements inwards from the boundary for 
     *             which all contained higher-order nodes will be free
     *             variables in the optimization.  Any vertex for further
     *             from the boundary will be slaved.  Specifying zero
     *             will result in all higher-order nodes being slaved
     *             except free nodes on the boundary.
     *\param max_boundary_domain_dimension  Specify the definition
     *             of "boundary".  If greater than or equal to 4, then
     *             the the set of all fixed vertices is assumed to be
     *             the boundary.  If less than four, then all vertices
     *             constrained to a domain with the specified number of
     *             fewer degrees of freedom (constrained to a geometric
     *             entity with an equal or smaller topological dimension)
     *             will be considered to be the boundary.
     */
    SlaveBoundaryVertices( unsigned depth,
                           unsigned max_boundary_domain_dimension = 4 );
    
    virtual double loop_over_mesh( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError& err );
   
    virtual std::string get_name() const;

    virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError& err );
    
    unsigned get_num_boundary_layers() const { return elemDepth; }
    bool boundary_is_fixed_vertices() const { return domainDoF >= 4; }
    bool boundary_is_mesh_domain() const { return domainDoF < 4; }
    unsigned boundary_mesh_domain_dimension() const { return domainDoF; }
    
   private:
    unsigned elemDepth;
    unsigned domainDoF;
};


} // namespace Mesquite

#endif
