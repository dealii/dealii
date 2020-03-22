/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   MeshTransform.hpp
  \brief  

  Class for performing an affine transformation on the mesh.

  \author Michael Brewer      
  \date   2004-11-06
*/

#ifndef Mesquite_MeshTransform_hpp 
#define Mesquite_MeshTransform_hpp


#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "Instruction.hpp"

namespace MESQUITE_NS
{

  /*! \class MeshTransform
    Perform an Affine transformation on Mesh vertex positions.
    Essentially define the new vertex position, v_new, from the original
    vertex position, v_old, s.t.
    v_new = (mMat * v_old) + mVec,
    where mMat is a constant matrix and mVec is a constant vector.
   */  
  class MESQUITE_EXPORT MeshTransform : public Instruction 
  {
  public:
    MeshTransform(bool skip_fixed = false) 
      : mMat(1,0,0,0,1,0,0,0,1), mVec(0.0), skipFixed(skip_fixed)
       {}
    MeshTransform(Matrix3D &in_mat, Vector3D &in_vec,
                                  bool skip_fixed = false)
      : mMat(in_mat), mVec(in_vec), skipFixed(skip_fixed)
      {}

      // virtual destructor ensures use of polymorphism during destruction
    virtual ~MeshTransform();
    
      //virtual functions from PatchDataUser...
      //!Loop over the mesh and perform the affine transformation
    virtual double loop_over_mesh( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError &err);
      //! Return the name of this PatchDataUser:  Mesh Transform
    virtual std::string get_name() const { return "Mesh Transform";}
    
    virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError& err );
    
    void add_translation( const Vector3D& offset );
    void add_rotation( const Vector3D& axis, double radians );
    void add_scale( double factor );
    void add_scale( const Vector3D& factors );
    
    bool skipping_fixed_vertices() const { return skipFixed; }
    void skip_fixed_vertices(bool yesno) { skipFixed = yesno; }
    
  private:
    Matrix3D mMat;//!Matrix for the affine transformation
    Vector3D mVec;//!Vector for the affine transformation
    bool skipFixed;
  };

  
} // namespace
#endif // Mesquite_MeshTransform_hpp
