/* $Id$ */
/*            Ralf Hartmann, University of Heidelberg, Dez 98               */

#include<fe/fe_lib.dg.h>



template<int dim>
FEDG_Q1<dim>::FEDG_Q1():
		FEQ1<dim>(1) {};


template<int dim>
FEDG_Q2<dim>::FEDG_Q2():
		FEQ2<dim>(1) {};


template<int dim>
FEDG_Q3<dim>::FEDG_Q3():
		FEQ3<dim>(1) {};


template<int dim>
FEDG_Q4<dim>::FEDG_Q4():
		FEQ4<dim>(1) {};



template <int dim>
void
FEDG_Q1<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
					  vector<Point<dim> >  &support_points) const {
  Assert ((support_points.size() == 0),
	  ExcWrongFieldDimension (support_points.size(),0));
};


template <int dim>
void
FEDG_Q2<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
						vector<Point<dim> >  &support_points) const {
  Assert ((support_points.size() == 0),
	  ExcWrongFieldDimension (support_points.size(),0));
};


template <int dim>
void
FEDG_Q3<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
					    vector<Point<dim> >  &support_points) const {
  Assert ((support_points.size() == 0),
	  ExcWrongFieldDimension (support_points.size(),0));
};


template <int dim>
void
FEDG_Q4<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
					      vector<Point<dim> >  &support_points) const {
  Assert ((support_points.size() == 0),
	  ExcWrongFieldDimension (support_points.size(),0));
};




// explicit instantiations
template class FEDG_Q1<deal_II_dimension>;
template class FEDG_Q2<deal_II_dimension>;
template class FEDG_Q3<deal_II_dimension>;
template class FEDG_Q4<deal_II_dimension>;
