
#include "SmartLaplacianSmoother.hpp"
#include "PatchData.hpp"
#include "MsqVertex.hpp"
#include "MsqMeshEntity.hpp"


namespace MESQUITE_NS {

size_t SmartLaplacianSmoother::num_inverted( PatchData& pd, MsqError& err )
{
  size_t result = 0;
  int inverted, junk;
  for (size_t i = 0; i < pd.num_elements(); ++i) {
    pd.element_by_index(i).check_element_orientation( pd, inverted, junk, err ); 
    MSQ_ERRZERO(err);
    if (inverted)
      ++result;
  }
  return result;
}


std::string SmartLaplacianSmoother::get_name() const
  { return "SmartLaplacianSmoother"; }

SmartLaplacianSmoother::~SmartLaplacianSmoother() 
{
}    

void SmartLaplacianSmoother::optimize_vertex_positions( PatchData &pd, 
                                                        MsqError &err )
{
  assert(pd.num_free_vertices() == 1);
  const size_t center_vtx_index = 0;
  const size_t init_inverted = num_inverted( pd, err ); MSQ_ERRRTN(err);
  
  adjVtxList.clear();
  pd.get_adjacent_vertex_indices( center_vtx_index, adjVtxList, err );
  MSQ_ERRRTN(err);
  
  if (adjVtxList.empty())
    return;
  
  const MsqVertex* verts = pd.get_vertex_array(err);
  const size_t n = adjVtxList.size();
  
  const Vector3D orig_pos = verts[center_vtx_index];
  Vector3D new_pos = verts[ adjVtxList[0] ];
  for (size_t i = 1; i < n; ++i)
    new_pos += verts[ adjVtxList[i] ];
  new_pos *= 1.0/n;
  pd.set_vertex_coordinates( new_pos, center_vtx_index, err );
  pd.snap_vertex_to_domain( center_vtx_index, err );  MSQ_ERRRTN(err);
  const size_t new_inverted = num_inverted( pd, err ); MSQ_ERRRTN(err);
  if (new_inverted > init_inverted)
    pd.set_vertex_coordinates( orig_pos, center_vtx_index, err );
}

}
