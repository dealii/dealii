/* $Id$ */
/*            Ralf Hartmann, University of Heidelberg, Dez 98               */

#include<fe/fe_lib.dg.h>



template<int dim>
FEDGLinear<dim>::FEDGLinear():
		FELinear<dim>(1) {};


template<int dim>
FEDGQuadraticSub<dim>::FEDGQuadraticSub():
		FEQuadraticSub<dim>(1) {};


template<int dim>
FEDGCubicSub<dim>::FEDGCubicSub():
		FECubicSub<dim>(1) {};


template<int dim>
FEDGQuarticSub<dim>::FEDGQuarticSub():
		FEQuarticSub<dim>(1) {};



template <int dim>
void
FEDGLinear<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
					  vector<Point<dim> >  &support_points) const {
  Assert ((support_points.size() == 0),
	  ExcWrongFieldDimension (support_points.size(),0));
};


template <int dim>
const FullMatrix<double> & 
FEDGLinear<dim>::restrict (const unsigned int child) const {
  Assert (false, ExcNotImplemented());
  return restriction[child];
};



template <int dim>
void
FEDGQuadraticSub<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
						vector<Point<dim> >  &support_points) const {
  Assert ((support_points.size() == 0),
	  ExcWrongFieldDimension (support_points.size(),0));
};


template <int dim>
const FullMatrix<double> & 
FEDGQuadraticSub<dim>::restrict (const unsigned int child) const {
  Assert (false, ExcNotImplemented());
  return restriction[child];
};



template <int dim>
void
FEDGCubicSub<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
					    vector<Point<dim> >  &support_points) const {
  Assert ((support_points.size() == 0),
	  ExcWrongFieldDimension (support_points.size(),0));
};


template <int dim>
const FullMatrix<double> & 
FEDGCubicSub<dim>::restrict (const unsigned int child) const {
  Assert (false, ExcNotImplemented());
  return restriction[child];
};




template <int dim>
void
FEDGQuarticSub<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
					      vector<Point<dim> >  &support_points) const {
  Assert ((support_points.size() == 0),
	  ExcWrongFieldDimension (support_points.size(),0));
};


template <int dim>
const FullMatrix<double> & 
FEDGQuarticSub<dim>::restrict (const unsigned int child) const {
  Assert (false, ExcNotImplemented());
  return restriction[child];
};




// explicit instantiations
template class FEDGLinear<deal_II_dimension>;
template class FEDGQuadraticSub<deal_II_dimension>;
template class FEDGCubicSub<deal_II_dimension>;
template class FEDGQuarticSub<deal_II_dimension>;
