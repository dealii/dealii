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


/** \file TargetWriter.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TargetWriter.hpp"
#include "TargetCalculator.hpp"
#include "WeightCalculator.hpp"
#include "MeshInterface.hpp"
#include "MsqMatrix.hpp"
#include "PatchData.hpp"
#include "PatchSet.hpp"
#include "MsqError.hpp"
#include "ElementPatches.hpp"
#include "ElemSampleQM.hpp"
#include <sstream>

namespace MESQUITE_NS {

TargetWriter::TargetWriter(  TargetCalculator* tc,
                             WeightCalculator* wc,
                             std::string target_base_name,
                             std::string weight_base_name )
  : targetCalc(tc), 
    weightCalc(wc), 
    targetName(target_base_name), 
    weightName(weight_base_name)
    {}

TargetWriter::~TargetWriter() {}

std::string TargetWriter::get_name() const
  { return "TargetWriter"; }

double TargetWriter::loop_over_mesh( MeshDomainAssoc* mesh_and_domain,
                                     const Settings* settings,
                                     MsqError& err )
{
  Mesh* mesh = mesh_and_domain->get_mesh();
  MeshDomain* domain = mesh_and_domain->get_domain();
 
  PatchData patch;
  patch.set_mesh( mesh );
  patch.set_domain( domain );
  if (settings)
    patch.attach_settings( settings );
  
  ElementPatches patch_set;
  patch_set.set_mesh( mesh );
  std::vector<PatchSet::PatchHandle> patches;
  std::vector<PatchSet::PatchHandle>::iterator p;
  std::vector<Mesh::VertexHandle> patch_verts;
  std::vector<Mesh::ElementHandle> patch_elems;
  
  patch_set.get_patch_handles( patches, err ); MSQ_ERRZERO(err);
  
  std::vector< MsqMatrix<3,3> > targets3d;
  std::vector< MsqMatrix<3,2> > targets2dorient;
  std::vector< MsqMatrix<2,2> > targets2d;
  std::vector< double > weights;
  std::vector< Sample > samples;
  for (p = patches.begin(); p != patches.end(); ++p)
  {
    patch_verts.clear();
    patch_elems.clear();
    patch_set.get_patch( *p, patch_elems, patch_verts, err ); MSQ_ERRZERO(err);
    patch.set_mesh_entities( patch_elems, patch_verts, err ); MSQ_ERRZERO(err);
    assert(patch.num_elements() == 1);
    
    MsqMeshEntity& elem = patch.element_by_index(0);
    EntityTopology type = elem.get_element_type();
    patch.get_samples( 0, samples, err ); MSQ_ERRZERO(err);
    if (samples.empty())
      continue;    
    
    if (targetCalc) {
      const unsigned dim = TopologyInfo::dimension(type);
      if (dim == 3) {
        targets3d.resize( samples.size() );
        for (unsigned i = 0; i < samples.size(); ++i) {
          targetCalc->get_3D_target( patch, 0, samples[i], targets3d[i], err ); MSQ_ERRZERO(err);

          if (DBL_EPSILON > det(targets3d[i])) {
            MSQ_SETERR(err)("Inverted 3D target", MsqError::INVALID_ARG);
            return 0.0;
          }
        }
        
        TagHandle tag = get_target_tag( 3, samples.size(), mesh,  err ); MSQ_ERRZERO(err);
        mesh->tag_set_element_data( tag, 1, 
                                    patch.get_element_handles_array(), 
                                    arrptr(targets3d), err ); MSQ_ERRZERO(err);
      }
      else if(targetCalc->have_surface_orient()) {
        targets2dorient.resize( samples.size() );
        for (unsigned i = 0; i < samples.size(); ++i) {
          targetCalc->get_surface_target( patch, 0, samples[i], targets2dorient[i], err ); MSQ_ERRZERO(err);

          MsqMatrix<3,1> cross = targets2dorient[i].column(0) * targets2dorient[i].column(1);
          if (DBL_EPSILON > (cross%cross)) {
            MSQ_SETERR(err)("Degenerate 2D target", MsqError::INVALID_ARG);
            return 0.0;
          }
        }
        
        TagHandle tag = get_target_tag( 2, samples.size(), mesh, err ); MSQ_ERRZERO(err);
        mesh->tag_set_element_data( tag, 1, 
                                    patch.get_element_handles_array(), 
                                    arrptr(targets2dorient), err ); MSQ_ERRZERO(err);
      }
      else {
        targets2d.resize( samples.size() );
        for (unsigned i = 0; i < samples.size(); ++i) {
          targetCalc->get_2D_target( patch, 0, samples[i], targets2d[i], err ); MSQ_ERRZERO(err);

          if (DBL_EPSILON > det(targets2d[i])) {
            MSQ_SETERR(err)("Degenerate/Inverted 2D target", MsqError::INVALID_ARG);
            return 0.0;
          }
        }
        
        TagHandle tag = get_target_tag( 2, samples.size(), mesh, err ); MSQ_ERRZERO(err);
        mesh->tag_set_element_data( tag, 1, 
                                    patch.get_element_handles_array(), 
                                    arrptr(targets2d), err ); MSQ_ERRZERO(err);
      }
    }
      
    if (weightCalc) {
      weights.resize( samples.size() );
      for (unsigned i = 0; i < samples.size(); ++i) {
        weights[i] = weightCalc->get_weight( patch, 0, samples[i], err ); MSQ_ERRZERO(err);
      }
      TagHandle tag = get_weight_tag( samples.size(), mesh, err ); MSQ_ERRZERO(err);
      mesh->tag_set_element_data( tag, 1, 
                                  patch.get_element_handles_array(),
                                  arrptr(weights), err ); MSQ_ERRZERO(err);
    }
  }

  return 0.0;
}

TagHandle TargetWriter::get_target_tag( unsigned dim, unsigned count, Mesh* mesh, MsqError& err )
{
  if (dim == 2 && !targetCalc->have_surface_orient())
    count *= 4;
  else
    count *= 3*dim; // num doubles
  if (targetTags.size() <= count)
    targetTags.resize( count+1, 0 );
  if (!targetTags[count]) {
    targetTags[count] = get_tag_handle( targetName, count, mesh, err );
    if (MSQ_CHKERR(err)) {
      targetTags[count] = 0;
      return 0;
    }
  }
  return targetTags[count];
}

TagHandle TargetWriter::get_weight_tag( unsigned count, Mesh* mesh, MsqError& err )
{
  if (weightTags.size() <= count)
    weightTags.resize( count+1, 0 );
  if (!weightTags[count]) {
    weightTags[count] = get_tag_handle( weightName, count, mesh, err );
    if (MSQ_CHKERR(err)) {
      weightTags[count] = 0;
      return 0;
    }
  }
  return weightTags[count];
}

TagHandle TargetWriter::get_tag_handle( const std::string& base_name,
                                        unsigned num_dbl, Mesh* mesh, MsqError& err )
{
  std::ostringstream sstr;
  sstr << base_name << num_dbl;
  
  TagHandle handle = mesh->tag_get( sstr.str().c_str(), err );
  if (!MSQ_CHKERR(err))
  {
    std::string temp_name;
    Mesh::TagType temp_type;
    unsigned temp_length;
    mesh->tag_properties( handle, temp_name, temp_type, temp_length, err );
    MSQ_ERRZERO(err);
    
    if (temp_type != Mesh::DOUBLE || temp_length != num_dbl)
    {
      MSQ_SETERR(err)( MsqError::TAG_ALREADY_EXISTS,
                      "Mismatched type or length for existing tag \"%s\"",
                       sstr.str().c_str() );
    }
  }
  else if (err.error_code() == MsqError::TAG_NOT_FOUND)
  {
    err.clear();
    handle = mesh->tag_create( sstr.str().c_str(), Mesh::DOUBLE, num_dbl, 0, err );
    MSQ_ERRZERO(err);
  }
  return handle;
}

void TargetWriter::initialize_queue( MeshDomainAssoc* ,
                                     const Settings* ,
                                     MsqError&  )
{
}
} // namespace Mesquite
