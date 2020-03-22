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


/** \file TargetCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TargetCalculator.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "ReferenceMesh.hpp"
#include "MappingFunction.hpp"
#include <assert.h>

namespace MESQUITE_NS {

double TargetCalculator::size( const MsqMatrix<3,3>& W )
{
  return Mesquite::cbrt( fabs(det(W)) );
}

double TargetCalculator::size( const MsqMatrix<3,2>& W )
{
  return sqrt(length(W.column(0) * W.column(1)));
}

double TargetCalculator::size( const MsqMatrix<2,2>& W )
{
  return sqrt( fabs(det(W)) );
}

MsqMatrix<3,3> TargetCalculator::skew( const MsqMatrix<3,3>& W )
{
  MsqVector<3> a1    = W.column(0) * (1.0/length(W.column(0)));
  MsqVector<3> a2    = W.column(1) * (1.0/length(W.column(1)));
  MsqVector<3> a3    = W.column(2) * (1.0/length(W.column(2)));
  MsqVector<3> a1xa2 = a1 * a2;
  
  double lenx  = length(a1xa2);
  double alpha = fabs(a1xa2 % a3);
  double coeff = Mesquite::cbrt(1/alpha);
  double dot   = a1xa2 % (a1 * a3);
  
  MsqMatrix<3,3> q;
  q(0,0) = coeff; q(0,1) = coeff * (a1 % a2); q(0,2) = coeff * (a1 % a3);
  q(1,0) = 0.0;   q(1,1) = coeff * lenx;      q(1,2) = coeff * dot   / lenx;
  q(2,0) = 0.0;   q(2,1) = 0.0;               q(2,2) = coeff * alpha / lenx;
  return q;
}


MsqMatrix<2,2> TargetCalculator::skew( const MsqMatrix<3,2>& W )
{
  MsqVector<3> alpha = W.column(0) * W.column(1);
  double a1_sqr      = W.column(0) % W.column(0);
  double a2_sqr      = W.column(1) % W.column(1);
  double dot         = W.column(0) % W.column(1);
  double a1a2        = sqrt(a1_sqr * a2_sqr);
  double coeff       = sqrt(a1a2/length(alpha));

  MsqMatrix<2,2> result;
  result(0,0) = coeff; result(0,1) = coeff * dot / a1a2;
  result(1,0) = 0.0  ; result(1,1) = 1/coeff;
  return result;
}


MsqMatrix<2,2> TargetCalculator::skew( const MsqMatrix<2,2>& W )
{
  double alpha   = fabs(det(W));
  double a1_sqr  = W.column(0) % W.column(0);
  double a2_sqr  = W.column(1) % W.column(1);
  double dot     = W.column(0) % W.column(1);
  double a1a2    = sqrt(a1_sqr * a2_sqr);
  double coeff   = sqrt(a1a2/alpha);

  MsqMatrix<2,2> result;
  result(0,0) = coeff; result(0,1) = coeff * dot / a1a2;
  result(1,0) = 0.0  ; result(1,1) = 1/coeff;
  return result;
}

MsqMatrix<3,3> TargetCalculator::aspect( const MsqMatrix<3,3>& W )
{
  double a1 = length(W.column(0));
  double a2 = length(W.column(1));
  double a3 = length(W.column(2));
  double coeff = 1.0/Mesquite::cbrt(a1*a2*a3);
  
  MsqMatrix<3,3> result(coeff);
  result(0,0) *= a1;
  result(1,1) *= a2;
  result(2,2) *= a3;
  return result;
}


MsqMatrix<2,2> TargetCalculator::aspect( const MsqMatrix<3,2>& W )
{
  double a1_sqr  = W.column(0) % W.column(0);
  double a2_sqr  = W.column(1) % W.column(1);
  double sqrt_rho = sqrt(sqrt(a1_sqr/a2_sqr));
  
  MsqMatrix<2,2> result;
  result(0,0) = sqrt_rho; result(0,1) = 0.0;
  result(1,0) = 0.0     ; result(1,1) = 1/sqrt_rho;
  return result;
}


MsqMatrix<2,2> TargetCalculator::aspect( const MsqMatrix<2,2>& W )
{
  double a1_sqr  = W.column(0) % W.column(0);
  double a2_sqr  = W.column(1) % W.column(1);
  double sqrt_rho = sqrt(sqrt(a1_sqr/a2_sqr));
  
  MsqMatrix<2,2> result;
  result(0,0) = sqrt_rho; result(0,1) = 0.0;
  result(1,0) = 0.0     ; result(1,1) = 1/sqrt_rho;
  return result;
}

MsqMatrix<3,3> TargetCalculator::shape( const MsqMatrix<3,3>& W )
{
  MsqVector<3> a1    = W.column(0);
  MsqVector<3> a2    = W.column(1);
  MsqVector<3> a3    = W.column(2);
  MsqVector<3> a1xa2 = a1 * a2;
  
  double len1  = length(a1);
  double lenx  = length(a1xa2);
  double alpha = fabs(a1xa2 % a3);
  double coeff = Mesquite::cbrt(1/alpha);
  double inv1  = 1.0/len1;
  double invx  = 1.0/lenx;
  
  MsqMatrix<3,3> q;
  q(0,0) = coeff*len1; q(0,1) = coeff*inv1*(a1 % a2); q(0,2) = coeff*inv1*(a1 % a3);
  q(1,0) = 0.0;        q(1,1) = coeff*inv1*lenx;      q(1,2) = coeff*invx*inv1*(a1xa2 % (a1 * a3));
  q(2,0) = 0.0;        q(2,1) = 0.0;                  q(2,2) = coeff * alpha * invx;
  return q;
}


MsqMatrix<2,2> TargetCalculator::shape( const MsqMatrix<3,2>& W )
{
  double len1 = length(W.column(0));
  double inv1 = 1.0/len1;
  double root_alpha = sqrt(length(W.column(0) * W.column(1)));
  double coeff   = 1.0/root_alpha;

  MsqMatrix<2,2> result;
  result(0,0) = coeff*len1; result(0,1) = coeff*inv1*(W.column(0) % W.column(1));
  result(1,0) = 0.0  ;      result(1,1) = root_alpha * inv1;
  return result;
}


MsqMatrix<2,2> TargetCalculator::shape( const MsqMatrix<2,2>& W )
{
  double len1 = length(W.column(0));
  double inv1 = 1.0/len1;
  double root_alpha = sqrt(fabs(det(W)));
  double coeff   = 1.0/root_alpha;

  MsqMatrix<2,2> result;
  result(0,0) = coeff*len1; result(0,1) = coeff*inv1*(W.column(0) % W.column(1));
  result(1,0) = 0.0       ; result(1,1) = root_alpha * inv1;
  return result;
}

bool TargetCalculator::factor_3D( const MsqMatrix<3,3>& A,
                                  double& Lambda,
                                  MsqMatrix<3,3>& V,
                                  MsqMatrix<3,3>& Q,
                                  MsqMatrix<3,3>& Delta,
                                  MsqError& err )
{
  MsqVector<3> a1xa2 = A.column(0) * A.column(1);
  double alpha = a1xa2 % A.column(2);
  Lambda = Mesquite::cbrt( fabs(alpha) );
  if (Lambda < DBL_EPSILON)
    return false;
  
  double la1_sqr = A.column(0) % A.column(0);
  double la1 = sqrt(la1_sqr);
  double la2 = length(A.column(1));
  double la3 = length(A.column(2));
  double lx = length(a1xa2);
  double a1dota2 = A.column(0) % A.column(1);
  
  double inv_la1 = 1.0/la1;
  double inv_lx  = 1.0/lx;
  V.set_column( 0, A.column(0) * inv_la1 );
  V.set_column( 1, (la1_sqr * A.column(1) - a1dota2 * A.column(0)) * inv_la1 * inv_lx );
  V.set_column( 2, (alpha / (fabs(alpha) * lx)) * a1xa2 );
  
  double inv_la2 = 1.0/la2;
  double inv_la3 = 1.0/la3;
  double len_prod_rt3 = Mesquite::cbrt( la1 * la2 * la3 );
  Q(0,0) = 1.0;
  Q(1,0) = Q(2,0) = 0.0;
  Q(0,1) = a1dota2 * inv_la1 * inv_la2;
  Q(1,1) = lx * inv_la1 * inv_la2;
  Q(2,1) = 0.0;
  Q(0,2) = (A.column(0) % A.column(2)) * inv_la1 * inv_la3;
  Q(1,2) = (a1xa2 % (A.column(0) * A.column(2))) * inv_lx * inv_la1 * inv_la3;
  Q(2,2) = fabs(alpha) * inv_lx * inv_la3;
  Q *= len_prod_rt3 / Lambda;
  
  double inv_prod_rt3 = 1.0/len_prod_rt3;;
  Delta(0,0) = la1*inv_prod_rt3;
  Delta(0,1) = 0.0;
  Delta(0,2) = 0.0;
  Delta(1,0) = 0.0;
  Delta(1,1) = la2*inv_prod_rt3;
  Delta(1,2) = 0.0;
  Delta(2,0) = 0.0;
  Delta(2,1) = 0.0;
  Delta(2,2) = la3*inv_prod_rt3;
  
  return true;
}

bool TargetCalculator::factor_surface( const MsqMatrix<3,2>& A,
                                  double& Lambda,
                                  MsqMatrix<3,2>& V,
                                  MsqMatrix<2,2>& Q,
                                  MsqMatrix<2,2>& Delta,
                                  MsqError& err )
{
  MsqVector<3> cross = A.column(0) * A.column(1);
  double alpha = length(cross);
  Lambda = sqrt(alpha);
  if (Lambda < DBL_EPSILON)
    return false;
  
  double la1_sqr = A.column(0) % A.column(0);
  double la1 = sqrt(la1_sqr);
  double la2 = length(A.column(1));
  double inv_la1 = 1.0/la1;
  double dot = A.column(0) % A.column(1);
  
  V.set_column( 0, A.column(0) * inv_la1 );
  V.set_column( 1, (la1_sqr * A.column(1) - dot * A.column(0)) / (la1*alpha) );
  
  double prod_rt2 = sqrt( la1 * la2 );
  Q(0,0) = prod_rt2 / Lambda; Q(0,1) = dot / (prod_rt2 * Lambda);
  Q(1,0) = 0.0; Q(1,1) = 1.0/Q(0,0);
  
  double inv_prod_rt2 = 1.0/prod_rt2;
  Delta(0,0) = la1*inv_prod_rt2; Delta(0,1) = 0.0;
  Delta(1,0) = 0.0;              Delta(1,1) = la2*inv_prod_rt2;
  
  return true;
}

bool TargetCalculator::factor_2D( const MsqMatrix<2,2>& A,
                                  double& Lambda,
                                  MsqMatrix<2,2>& V,
                                  MsqMatrix<2,2>& Q,
                                  MsqMatrix<2,2>& Delta,
                                  MsqError& err )
{
  double alpha = det(A);
  Lambda = sqrt(fabs(alpha));
  if (Lambda < DBL_EPSILON)
    return false;
  
  double la1_sqr = A.column(0) % A.column(0);
  double la1 = sqrt(la1_sqr);
  double la2 = length(A.column(1));
  double inv_la1 = 1.0/la1;
  double dot = A.column(0) % A.column(1);
  
  V.set_column( 0, A.column(0) * inv_la1 );
  V.set_column( 1, (la1_sqr * A.column(1) - dot * A.column(0)) / (la1*alpha) );
  
  double prod_rt2 = sqrt( la1 * la2 );
  Q(0,0) = prod_rt2 / Lambda; Q(0,1) = dot / (prod_rt2 * Lambda);
  Q(1,0) = 0.0; Q(1,1) = 1.0/Q(0,0);
  
  double inv_prod_rt2 = 1.0/prod_rt2;
  Delta(0,0) = la1*inv_prod_rt2; Delta(0,1) = 0.0;
  Delta(1,0) = 0.0;              Delta(1,1) = la2*inv_prod_rt2;
  
  return true;
}

MsqMatrix<3,3> TargetCalculator::new_orientation_3D( const MsqVector<3>& b1,
                                                     const MsqVector<3>& b2 )
{
  double lb1_sqr = b1 % b1;
  MsqVector<3> cross = b1 * b2;
  double lb1 = sqrt(lb1_sqr);
  double inv_lb1 = 1.0/lb1;
  double inv_lx = 1.0/length(cross);
  MsqMatrix<3,3> V;
  V.set_column( 0, inv_lb1 * b1 );
  V.set_column( 1, (lb1_sqr * b2 - (b1 % b2) * lb1) * inv_lb1 * inv_lx );
  V.set_column( 2, cross * inv_lx );
  return V;
}
                                                     
MsqMatrix<3,2> TargetCalculator::new_orientation_2D( const MsqVector<3>& b1,
                                                     const MsqVector<3>& b2 )
{
  double lb1_sqr = b1 % b1;
  double inv_lb1 = 1.0/sqrt(lb1_sqr);
  MsqMatrix<3,2> V;
  V.set_column(0, b1 * inv_lb1);
  V.set_column(1, (lb1_sqr * b2 - (b1 % b2) * b1) * (inv_lb1 / length(b1 * b2)));
  return V;
}

/** If, for the specified element type, the skew is constant for
 *  an ideal element and the aspect is identity everywhere within
 *  the element, pass back the constant skew/shape term and return 
 *  true.  Otherwise return false.
 */
static inline bool ideal_constant_skew_I_3D( EntityTopology element_type,
                                      MsqMatrix<3,3>& q )
{
  switch (element_type) {
    case HEXAHEDRON:
      q = MsqMatrix<3,3>(1.0); // Identity
      return false;
    case TETRAHEDRON:
      // [ x,   x/2, x/2 ] x^6 = 2
      // [ 0,   y,   y/3 ] y^2 = 3/4 x^2
      // [ 0,   0,   z   ] z^2 = 2/3 x^2
      q(0,0) = 1.122462048309373; 
      q(0,1) = q(0,2) = 0.56123102415468651;
      q(1,0) = q(2,0) = q(2,1) = 0.0;
      q(1,1) = 0.97208064861983279;
      q(1,2) = 0.32402688287327758;
      q(2,2) = 0.91648642466573493;
      return true;
    case PRISM:
      //            [ 1   0   0 ]
      //  a^(-1/3)  [ 0   1  1/2]
      //            [ 0   0   a ]
      //
      // a = sqrt(3)/2
      //
      q(0,0) = q(1,1) = 1.0491150634216482;
      q(0,1) = q(0,2) = q(1,0) = 0.0;
      q(1,2) = 0.52455753171082409;
      q(2,0) = q(2,1) = 0.0;
      q(2,2) = 0.90856029641606972;
      return true;
    default:
      return false;
  }
}

/** If, for the specified element type, the skew is constant for
 *  an ideal element and the aspect is identity everywhere within
 *  the element, pass back the constant skew/shape term and return 
 *  true.  Otherwise return false.
 */
static inline bool ideal_constant_skew_I_2D( EntityTopology element_type,
                                      MsqMatrix<2,2>& q )
{
  switch (element_type) {
    case QUADRILATERAL:
      q = MsqMatrix<2,2>(1.0); // Identity
      return true;
    case TRIANGLE:
      // [ x, x/2 ]  x = pow(4/3, 0.25)
      // [ 0, y   ]  y = 1/x
      q(0,0) = 1.074569931823542; q(0,1) = 0.537284965911771;
      q(1,0) = 0.0;               q(1,1) = 0.93060485910209956;
      return true;
    default:
      return false;
  }
}

void TargetCalculator::ideal_skew_3D( EntityTopology element_type,
                                      Sample s,
                                      const PatchData& pd,
                                      MsqMatrix<3,3>& q,
                                      MsqError& err ) 
{
  if (!ideal_constant_skew_I_3D(element_type,q)) {
    const MappingFunction3D* map = pd.get_mapping_function_3D( element_type );
    if (!map) {
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
      return;
    }
    map->ideal( s, q, err );
    MSQ_ERRRTN(err);
    q = TargetCalculator::skew(q);
  }
}

void TargetCalculator::ideal_skew_2D( EntityTopology element_type,
                                      Sample s,
                                      const PatchData& pd,
                                      MsqMatrix<2,2>& q,
                                      MsqError& err ) 
{
   if (!ideal_constant_skew_I_2D(element_type,q)) {
    const MappingFunction2D* map = pd.get_mapping_function_2D( element_type );
    if (!map) {
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
      return;
    }
    MsqMatrix<3,2> J;
    map->ideal( s, J, err );
    MSQ_ERRRTN(err);
    q = TargetCalculator::skew(J);
  }
}

void TargetCalculator::ideal_shape_3D( EntityTopology element_type,
                                      Sample s,
                                      const PatchData& pd,
                                      MsqMatrix<3,3>& q,
                                      MsqError& err ) 
{
  if (!ideal_constant_skew_I_3D(element_type,q)) {
    const MappingFunction3D* map = pd.get_mapping_function_3D( element_type );
    if (!map) {
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
      return;
    }
    map->ideal( s, q, err );
    MSQ_ERRRTN(err);
    q = TargetCalculator::shape(q);
  }
}

void TargetCalculator::ideal_shape_2D( EntityTopology element_type,
                                      Sample s,
                                      const PatchData& pd,
                                      MsqMatrix<2,2>& q,
                                      MsqError& err ) 
{
   if (!ideal_constant_skew_I_2D(element_type,q)) {
    const MappingFunction2D* map = pd.get_mapping_function_2D( element_type );
    if (!map) {
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
      return;
    }
    MsqMatrix<3,2> J;
    map->ideal( s, J, err );
    MSQ_ERRRTN(err);
    q = TargetCalculator::shape(J);
  }
}

MsqMatrix<3,3> TargetCalculator::new_aspect_3D( const MsqVector<3>& r )
{
  MsqMatrix<3,3> W(0.0);
  W(0,0) = Mesquite::cbrt( r[0]*r[0]/(r[1]*r[2]) );
  W(1,1) = Mesquite::cbrt( r[1]*r[1]/(r[0]*r[2]) );
  W(2,2) = Mesquite::cbrt( r[2]*r[2]/(r[1]*r[0]) );
  return W;
}

MsqMatrix<2,2> TargetCalculator::new_aspect_2D( const MsqVector<2>& r )
{
  return new_aspect_2D( r[0] / r[1] );
}

MsqMatrix<2,2> TargetCalculator::new_aspect_2D( double rho )
{
  MsqMatrix<2,2> W( sqrt(rho) );
  W(1,1) = 1.0/W(0,0);
  return W;
}

static NodeSet get_nodeset( EntityTopology type, int num_nodes, MsqError& err )
{
  bool midedge, midface, midvol;
  TopologyInfo::higher_order( type, num_nodes, midedge, midface, midvol, err );
  if (MSQ_CHKERR(err)) return NodeSet();
  
  NodeSet bits;
  bits.set_all_corner_nodes(type);
  if (midedge) 
    bits.set_all_mid_edge_nodes(type);
  if (midface)
    bits.set_all_mid_face_nodes(type);
  if (TopologyInfo::dimension(type) == 3 && midvol)
    bits.set_mid_region_node();
 
  return bits;
}

void TargetCalculator::jacobian_3D( PatchData& pd, 
                                    EntityTopology type,
                                    int num_nodes,
                                    Sample location,
                                    const Vector3D* coords,
                                    MsqMatrix<3,3>& J,
                                    MsqError& err )
{
    // Get element properties
  NodeSet bits = get_nodeset( type, num_nodes, err ); MSQ_ERRRTN(err);
  const MappingFunction3D* mf = pd.get_mapping_function_3D( type );
  if (!mf) {
    MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
    return;
  }

    // Get mapping function derivatives
  const int MAX_NODES = 27;
  assert(num_nodes <= MAX_NODES);
  size_t indices[MAX_NODES], n;
  MsqVector<3> derivs[MAX_NODES];
  mf->derivatives( location, bits, indices, derivs, n, err ); MSQ_ERRRTN(err);
  
    // calculate Jacobian
  assert(sizeof(Vector3D) == sizeof(MsqVector<3>));
  const MsqVector<3>* verts = reinterpret_cast<const MsqVector<3>*>(coords);
  assert(n > 0);
  J = outer( verts[indices[0]], derivs[0]  );
  for (size_t i = 1; i < n; ++i) 
    J += outer( verts[indices[i]], derivs[i] );
}

void TargetCalculator::jacobian_2D( PatchData& pd,
                                    EntityTopology type,
                                    int num_nodes,
                                    Sample location,
                                    const Vector3D* coords,
                                    MsqMatrix<3,2>& J,
                                    MsqError& err )
{
    // Get element properties
  NodeSet bits = get_nodeset( type, num_nodes, err ); MSQ_ERRRTN(err);
  const MappingFunction2D* mf = pd.get_mapping_function_2D( type );
  if (!mf) {
    MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
    return;
  }

    // Get mapping function derivatives
  const int MAX_NODES = 9;
  assert(num_nodes <= MAX_NODES);
  size_t indices[MAX_NODES], n;
  MsqVector<2> derivs[MAX_NODES];
  mf->derivatives( location, bits, indices, derivs, n, err ); MSQ_ERRRTN(err);
  
    // calculate Jacobian
  assert(sizeof(Vector3D) == sizeof(MsqVector<3>));
  const MsqVector<3>* verts = reinterpret_cast<const MsqVector<3>*>(coords);
  assert(n > 0);
  J = outer( verts[indices[0]], derivs[0] );
  for (size_t i = 1; i < n; ++i) 
    J += outer( verts[indices[i]], derivs[i] );
}


void TargetCalculator::get_refmesh_Jacobian_3D( 
                              ReferenceMeshInterface* ref_mesh,
                              PatchData& pd,
                              size_t element,
                              Sample sample,
                              MsqMatrix<3,3>& W_out,
                              MsqError& err )
{
    // get element
  MsqMeshEntity& elem = pd.element_by_index( element );
  const EntityTopology type = elem.get_element_type();
  const unsigned n = elem.node_count();

  const unsigned MAX_NODES = 27;
  assert(n <= MAX_NODES);
  
    // get vertices
  Mesh::VertexHandle elem_verts[MAX_NODES];
  const std::size_t* vtx_idx = elem.get_vertex_index_array();
  const Mesh::VertexHandle* vtx_hdl = pd.get_vertex_handles_array();
  for (unsigned i = 0; i < n; ++i)
    elem_verts[i] = vtx_hdl[vtx_idx[i]];
  
    // get vertex coordinates
  Vector3D vert_coords[MAX_NODES];
  ref_mesh->get_reference_vertex_coordinates( elem_verts, n, vert_coords, err );
  MSQ_ERRRTN(err);
  
    // calculate Jacobian
  jacobian_3D( pd, type, n, sample, vert_coords, W_out, err );
  MSQ_ERRRTN(err);
}

void TargetCalculator::get_refmesh_Jacobian_2D( 
                              ReferenceMeshInterface* ref_mesh,
                              PatchData& pd,
                              size_t element,
                              Sample sample,
                              MsqMatrix<3,2>& W_out,
                              MsqError& err )
{
    // get element
  MsqMeshEntity& elem = pd.element_by_index( element );
  const EntityTopology type = elem.get_element_type();
  const unsigned n = elem.node_count();

  const unsigned MAX_NODES = 9;
  assert(n <= MAX_NODES);
  
    // get vertices
  Mesh::VertexHandle elem_verts[MAX_NODES];
  const std::size_t* vtx_idx = elem.get_vertex_index_array();
  const Mesh::VertexHandle* vtx_hdl = pd.get_vertex_handles_array();
  for (unsigned i = 0; i < n; ++i)
    elem_verts[i] = vtx_hdl[vtx_idx[i]];
  
    // get vertex coordinates
  Vector3D vert_coords[MAX_NODES];
  ref_mesh->get_reference_vertex_coordinates( elem_verts, n, vert_coords, err );
  MSQ_ERRRTN(err);
  
    // calculate Jacobian
  jacobian_2D( pd, type, n, sample, vert_coords, W_out, err );
  MSQ_ERRRTN(err);
}

TargetCalculator::~TargetCalculator()
{}

void TargetCalculator::initialize_queue( MeshDomainAssoc* ,
                                         const Settings* ,
                                         MsqError&  )
{}


} // namespace MESQUITE_NS
