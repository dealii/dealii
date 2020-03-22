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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TargetMetricUtil.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TargetMetricUtil.hpp"
#include "MsqMatrix.hpp"
#include "PatchData.hpp"
#include "ElementQM.hpp"
#include "ElemSampleQM.hpp"

namespace MESQUITE_NS {

void surface_to_2d( const MsqMatrix<3,2>& A,
                    const MsqMatrix<3,2>& W,
                    MsqMatrix<2,2>& W_22,
                    MsqMatrix<3,2>& RZ )
{
  MsqMatrix<3,1> W1 = W.column(0);
  MsqMatrix<3,1> W2 = W.column(1);
  MsqMatrix<3,1> nw = W1 * W2;
  nw *= 1.0/length(nw);

  MsqMatrix<3,1> z[2];
  z[0] = W1 * (1.0 / length( W1 ));
  z[1] = nw * z[0];
  MsqMatrix<3,2> Z(z);

  MsqMatrix<3,1> np = A.column(0) * A.column(1);
  np *= 1.0 / length(np);
  double dot = np % nw;
  MsqMatrix<3,1> nr = (dot >= 0.0) ? nw : -nw;
  MsqMatrix<3,1> v = nr * np;
  double vlen = length(v);
  if (vlen < DBL_EPSILON) {
    RZ = Z; // R = I
  }
  else {
    v *= 1.0 / vlen;
    MsqMatrix<3,1> r1[3] = { v, np, v * np };
    MsqMatrix<3,1> r2[3] = { v, nr, v * nr };
    MsqMatrix<3,3> R1( r1 ), R2( r2 );
    RZ = R1 * transpose(R2) * Z;
  }
  
  W_22 = transpose(Z) * W;
}
/*
void surface_to_2d( const MsqMatrix<3,2>& App,
                    const MsqMatrix<3,2>& Wp,
                    MsqMatrix<2,2>& A,
                    MsqMatrix<2,2>& W )
{
  MsqMatrix<3,1> Wp1 = Wp.column(0);
  MsqMatrix<3,1> Wp2 = Wp.column(1);
  MsqMatrix<3,1> nwp = Wp1 * Wp2;
  nwp *= 1.0/length(nwp);

  MsqMatrix<3,1> z[2];
  z[0] = Wp1 * (1.0 / length( Wp1 ));
  z[1] = nwp * z[0];
  MsqMatrix<3,2> Z(z);
  W = transpose(Z) * Wp;

  MsqMatrix<3,1> npp = App.column(0) * App.column(1);
  npp *= 1.0 / length(npp);
  double dot = npp % nwp;
  MsqMatrix<3,1> nr = (dot >= 0.0) ? nwp : -nwp;
  MsqMatrix<3,1> v = nr * npp;
  double vlen = length(v);
  if (vlen > DBL_EPSILON) {
    v *= 1.0 / vlen;
    MsqMatrix<3,1> r1[3] = { v, npp, v * npp }, r2[3] = { v, nr, v * nr };
    MsqMatrix<3,3> R1( r1 ), R2( r2 );
    MsqMatrix<3,3> RT = R2 * transpose(R1);
    MsqMatrix<3,2> Ap = RT * App;
    A = transpose(Z) * Ap;
  }
  else {
    A = transpose(Z) * App;
  }
}
*/

static inline void append_elem_samples( PatchData& pd,
                                        size_t e,
                                        std::vector<size_t>& handles )
{
  NodeSet samples = pd.get_samples( e );
  EntityTopology type = pd.element_by_index( e ).get_element_type();
  size_t curr_size = handles.size();
  handles.resize( curr_size + samples.num_nodes() );
  std::vector<size_t>::iterator i = handles.begin() + curr_size;
  for (unsigned j = 0; j < TopologyInfo::corners(type); ++j)
    if (samples.corner_node(j))
      *(i++) = ElemSampleQM::handle( Sample(0,j), e );
  for (unsigned j = 0; j < TopologyInfo::edges(type); ++j)
    if (samples.mid_edge_node(j))
      *(i++) = ElemSampleQM::handle( Sample(1,j), e );
  for (unsigned j = 0; j < TopologyInfo::faces(type); ++j)
    if (samples.mid_face_node(j))
      *(i++) = ElemSampleQM::handle( Sample(2,j), e );
  if (TopologyInfo::dimension(type) == 3 && samples.mid_region_node())
    *(i++) = ElemSampleQM::handle( Sample(3,0), e );
}  
                                        

void get_sample_pt_evaluations( PatchData& pd,
                                std::vector<size_t>& handles,
                                bool free,
                                MsqError& err )
{
  handles.clear();
  std::vector<size_t> elems;
  ElementQM::get_element_evaluations( pd, elems, free, err ); MSQ_ERRRTN(err);
  for (std::vector<size_t>::iterator i = elems.begin(); i != elems.end(); ++i)
    append_elem_samples( pd, *i, handles );
}
                   
void get_elem_sample_points( PatchData& pd,
                             size_t elem,
                             std::vector<size_t>& handles,
                             MsqError& err )
{
  handles.clear();
  append_elem_samples( pd, elem, handles );
}

} // namespace Mesquite
