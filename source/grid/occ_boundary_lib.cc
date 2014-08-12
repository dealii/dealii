#include <deal.II/grid/occ_boundary_lib.h>

#ifdef DEAL_II_WITH_OPENCASCADE

DEAL_II_NAMESPACE_OPEN

namespace OpenCASCADE
{
  template <int dim, int spacedim>
  NormalProjectionBoundary<dim,spacedim>::NormalProjectionBoundary(const TopoDS_Shape sh, 
								   const double tolerance) :
    sh(sh),
    tolerance(tolerance) 
  {
    Assert(spacedim == 3, ExcImpossibleInDim(spacedim));
  }
  
  
  template <int dim, int spacedim>
  Point<spacedim>  NormalProjectionBoundary<dim,spacedim>::
  project_to_manifold (const std::vector<Point<spacedim> > &surrounding_points,
		       const Point<spacedim> &candidate) const {
    TopoDS_Shape out_shape;
    double u=0, v=0;
    for(unsigned int i=0; i<surrounding_points.size(); ++i)
      Assert(closest_point(sh, surrounding_points[i], out_shape, u, v)
	     .distance(surrounding_points[i]) < (surrounding_points[i].norm()>0 ? 
						 tolerance*surrounding_points[i].norm() :
						 tolerance),
	     ExcPointNotOnManifold(surrounding_points[i]));
    
    return closest_point(sh, candidate, out_shape, u, v);
  }
  
  // Explicit instantiations
#include "occ_boundary_lib.inst"  
  
} // end namespace OpenCASCADE

DEAL_II_NAMESPACE_CLOSE

#endif
