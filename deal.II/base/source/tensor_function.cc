// $Id$


#include <base/tensor_function.h>
#include <vector>
#include <base/tensor.h>
#include <cmath>
#include <lac/vector.h>



//////////////////////////////////////////////////////////////////////
// TensorFunction
//////////////////////////////////////////////////////////////////////

template <int rank, int dim>
TensorFunction<rank, dim>::TensorFunction (const double initial_time)
		:
		FunctionTime (initial_time)
{};



template <int rank, int dim>
TensorFunction<rank, dim>::~TensorFunction ()
{};



template <int rank, int dim>
Tensor<rank,dim>
TensorFunction<rank, dim>::value (const Point<dim> &) const
{
  Assert (false, ExcPureFunctionCalled());
  return Tensor<rank,dim>();
};



template <int rank, int dim>
void
TensorFunction<rank, dim>::value_list (const vector<Point<dim> > &points,
				       vector<Tensor<rank,dim> > &values) const
{
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    values[i]  = this->value (points[i]);
};



template <int rank, int dim>
Tensor<rank+1,dim>
TensorFunction<rank, dim>::gradient (const Point<dim> &) const
{
  Assert (false, ExcPureFunctionCalled());
  return Tensor<rank+1,dim>();
};



template <int rank, int dim>
void
TensorFunction<rank, dim>::gradient_list (const vector<Point<dim> > &points,
				    vector<Tensor<rank+1,dim> > &gradients) const
{
  Assert (gradients.size() == points.size(),
	  ExcVectorHasWrongSize(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i] = gradient(points[i]);
};



template class TensorFunction<1,1>;
template class TensorFunction<2,1>;
template class TensorFunction<3,1>;
template class TensorFunction<4,1>;
template class TensorFunction<1,2>;
template class TensorFunction<2,2>;
template class TensorFunction<3,2>;
template class TensorFunction<4,2>;
template class TensorFunction<1,3>;
template class TensorFunction<2,3>;
template class TensorFunction<3,3>;
template class TensorFunction<4,3>;
