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
