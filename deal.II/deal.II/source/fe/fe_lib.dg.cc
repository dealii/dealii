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



#if deal_II_dimension == 1

template <>
const FiniteElementData<1>
FEDGLinear<1>::get_fe_data () {
				   // no dofs at the vertices, all in the
				   // interior of the line
  return FiniteElementData<1> (0,
			       FELinear<1>::get_fe_data().total_dofs,
			       FELinear<1>::get_fe_data().n_transform_functions);
};


template <>
const FiniteElementData<1>
FEDGQuadraticSub<1>::get_fe_data () {
				   // no dofs at the vertices, all in the
				   // interior of the line
  return FiniteElementData<1> (0,
			       FEQuadraticSub<1>::get_fe_data().total_dofs,
			       FEQuadraticSub<1>::get_fe_data().n_transform_functions);
};


template <>
const FiniteElementData<1>
FEDGCubicSub<1>::get_fe_data () {
				   // no dofs at the vertices, all in the
				   // interior of the line
  return FiniteElementData<1> (0,
			       FECubicSub<1>::get_fe_data().total_dofs,
			       FECubicSub<1>::get_fe_data().n_transform_functions);
};


template <>
const FiniteElementData<1>
FEDGQuarticSub<1>::get_fe_data () {
				   // no dofs at the vertices, all in the
				   // interior of the line
  return FiniteElementData<1> (0,
			       FEQuarticSub<1>::get_fe_data().total_dofs,
			       FEQuarticSub<1>::get_fe_data().n_transform_functions);
};

#endif



#if deal_II_dimension == 2

template <>
const FiniteElementData<2>
FEDGLinear<2>::get_fe_data () {
				   // no dofs at the vertices or lines, all in the
				   // interior of the line
  return FiniteElementData<2> (0, 0,
			       FELinear<2>::get_fe_data().total_dofs,
			       FELinear<2>::get_fe_data().n_transform_functions);
};


template <>
const FiniteElementData<2>
FEDGQuadraticSub<2>::get_fe_data () {
				   // no dofs at the vertices or lines, all in the
				   // interior of the line
  return FiniteElementData<2> (0, 0,
			       FEQuadraticSub<2>::get_fe_data().total_dofs,
			       FEQuadraticSub<2>::get_fe_data().n_transform_functions);
};


template <>
const FiniteElementData<2>
FEDGCubicSub<2>::get_fe_data () {
				   // no dofs at the vertices or lines, all in the
				   // interior of the line
  return FiniteElementData<2> (0, 0,
			       FECubicSub<2>::get_fe_data().total_dofs,
			       FECubicSub<2>::get_fe_data().n_transform_functions);
};


template <>
const FiniteElementData<2>
FEDGQuarticSub<2>::get_fe_data () {
				   // no dofs at the vertices or lines, all in the
				   // interior of the line
  return FiniteElementData<2> (0, 0,
			       FEQuarticSub<2>::get_fe_data().total_dofs,
			       FEQuarticSub<2>::get_fe_data().n_transform_functions);
};

#endif



template <int dim>
const dFMatrix & 
FEDGLinear<dim>::restrict (const unsigned int child) const {
  Assert (false, ExcNotImplemented());
  return restriction[child];
};


template <int dim>
const dFMatrix & 
FEDGQuadraticSub<dim>::restrict (const unsigned int child) const {
  Assert (false, ExcNotImplemented());
  return restriction[child];
};



template <int dim>
const dFMatrix & 
FEDGCubicSub<dim>::restrict (const unsigned int child) const {
  Assert (false, ExcNotImplemented());
  return restriction[child];
};



template <int dim>
const dFMatrix & 
FEDGQuarticSub<dim>::restrict (const unsigned int child) const {
  Assert (false, ExcNotImplemented());
  return restriction[child];
};




// explicit instantiations
template class FEDGLinear<deal_II_dimension>;
template class FEDGQuadraticSub<deal_II_dimension>;
template class FEDGCubicSub<deal_II_dimension>;
template class FEDGQuarticSub<deal_II_dimension>;
