/*            Ralf Hartmann, University of Heidelberg, Dez 98               */


#include <fe/fe_lib.dg.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>
#include <algorithm>





template <int dim>
FEDG_Q0<dim>::FEDG_Q0 () :
		FEQ1Mapping<dim> (0, 
				  (dim==1 ? 1 : 0),
				  (dim==2 ? 1 : 0),
				  (dim==3 ? 1 : 0),
				  1,
				  true)
{
  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i)
    { 
      restriction[i](0,0) = 1./GeometryInfo<dim>::children_per_cell;
      prolongation[i](0,0) = 1.0;
    }
};



#if deal_II_dimension == 1


template <>
void
FEDG_Q0<1>::get_face_support_points (const DoFHandler<1>::face_iterator &,
				     vector<Point<1> >  &) const {
  Assert (false, ExcInternalError());
};

#endif


template <int dim>
inline
double
FEDG_Q0<dim>::shape_value (const unsigned int i,
				const Point<dim>&) const
{
  Assert((i<total_dofs), ExcIndexRange(i, 0, total_dofs));
  return 1.;
};



template <int dim>
inline
Tensor<1,dim>
FEDG_Q0<dim>::shape_grad (const unsigned int i,
			       const Point<dim>&) const
{
  Assert((i<total_dofs), ExcIndexRange(i, 0, total_dofs));
  return Tensor<1,dim> ();
};



template <int dim>
inline
Tensor<2,dim>
FEDG_Q0<dim>::shape_grad_grad (const unsigned int i,
				  const Point<dim> &) const
{
  Assert((i<total_dofs), ExcIndexRange(i, 0, total_dofs));

  return Tensor<2,dim>();
};



template <int dim>
void FEDG_Q0<dim>::get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					       FullMatrix<double> &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

  local_mass_matrix(0,0) = cell->measure();
};



template <int dim>
void
FEDG_Q0<dim>::get_unit_support_points (vector<Point<dim> > &unit_points) const {
  Assert (unit_points.size() == total_dofs,
	  ExcWrongFieldDimension (unit_points.size(), total_dofs));
  for (unsigned int d=0; d<dim; ++d)
    unit_points[0](d) = 0.5;
};
  

template <int dim>
void
FEDG_Q0<dim>::get_support_points (const typename DoFHandler<dim>::cell_iterator &cell,
				       vector<Point<dim> >  &support_points) const {
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));
  
  support_points[0] = cell->center();
};


template <int dim>
void
FEDG_Q0<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
					    vector<Point<dim> >  &support_points) const {
  Assert ((support_points.size() == 0),
	  ExcWrongFieldDimension (support_points.size(),0));
};



// explicit instantiations

template class FEDG_Q0<deal_II_dimension>;
