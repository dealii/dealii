//----------------------------  tensor_function.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_function.cc  ---------------------------


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


// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif


template <int rank, int dim>
void
TensorFunction<rank, dim>::value_list (const typename std::vector<Point<dim> > &points,
				       typename std::vector<Tensor<rank,dim> > &values) const
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch(values.size(), points.size()));

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
TensorFunction<rank, dim>::gradient_list (const typename std::vector<Point<dim> >   &points,
					  typename std::vector<Tensor<rank+1,dim> > &gradients) const
{
  Assert (gradients.size() == points.size(),
	  ExcDimensionMismatch(gradients.size(), points.size()));

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
