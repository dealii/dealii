
#ifndef Mesquite_SmartLaplacianSmoother_hpp 
#define Mesquite_SmartLaplacianSmoother_hpp

#include "Mesquite.hpp"
#include "RelaxationSmoother.hpp"

#include <vector>

namespace MESQUITE_NS
{
  /*\brief Do laplacian smooth, but don't invert elements.
   */  
  class SmartLaplacianSmoother : public RelaxationSmoother 
  {
  public:
    /**
     *\param OF ObjectiveFunction used by some termination criteria
     */
	MESQUITE_EXPORT
    SmartLaplacianSmoother( ObjectiveFunction* OF = NULL ) 
      : RelaxationSmoother(OF) {}
    
	MESQUITE_EXPORT
    ~SmartLaplacianSmoother();

	MESQUITE_EXPORT
    virtual std::string get_name() const;
    
	MESQUITE_EXPORT
    static size_t num_inverted( PatchData& pd, MsqError& err );
    
  protected:
	MESQUITE_EXPORT
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);

  private:
    std::vector<size_t> adjVtxList;    
  };

  

  
}

#endif
