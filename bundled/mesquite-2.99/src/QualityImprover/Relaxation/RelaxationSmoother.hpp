#ifndef RELAXATION_SMOOTHER_HPP
#define RELAXATION_SMOOTHER_HPP

#include "Mesquite.hpp"
#include "VertexPatches.hpp"
#include "VertexMover.hpp"

namespace MESQUITE_NS
{
  /**\brief Base class for LaPlacian and other relaxation smoothers */
  class RelaxationSmoother : public VertexMover
  {
  public:
    /**
     * \param OF  For many relaxation solvers (e.g. Laplacian)
     * this ObjectiveFunction is used only to evaluate user-specified
     * termination criteria that require an objective function.  
     */
	MESQUITE_EXPORT
    RelaxationSmoother( ObjectiveFunction* OF = NULL ) : VertexMover(OF) {}
    
	MESQUITE_EXPORT
    virtual ~RelaxationSmoother();

	MESQUITE_EXPORT
    PatchSet* get_patch_set() { return &patchSet; }
    
  protected:
	MESQUITE_EXPORT
    virtual void initialize(PatchData &pd, MsqError &err);

	MESQUITE_EXPORT
    virtual void optimize_vertex_positions( PatchData &pd,
                                            MsqError &err) = 0;
					 
	MESQUITE_EXPORT
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);

	MESQUITE_EXPORT
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);

	MESQUITE_EXPORT
    virtual void cleanup();
    
  private:
    VertexPatches patchSet;
  };
  
}

#endif
