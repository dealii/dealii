// $Id$


#include <base/tensorfunction.h>
#include <vector>
#include <base/tensor.h>
#include <cmath>

template <int dim>
VectorFunction<dim>::VectorFunction(unsigned n_components, const double initial_time)
		:
		FunctionTime(initial_time),
		n_components(n_components)
{}


template <int dim>
VectorFunction<dim>::~VectorFunction()
{}

/*
template <int dim> double
VectorFunction<dim>::operator () (const Point<dim> &, unsigned) const

{
  Assert (false, ExcPureFunctionCalled());
  return 0.;
}
*/

template <int dim> void
VectorFunction<dim>::value (const Point<dim>  &, vector<double> &) const
{
  Assert (false, ExcPureFunctionCalled());
}


template <int dim> void
VectorFunction<dim>::value_list (const vector<Point<dim> > &,
				 vector<vector<double> > &) const
{
  Assert (false, ExcPureFunctionCalled());
}


template <int dim> void
VectorFunction<dim>::gradient_list (const vector<Point<dim> > &,
				    vector<vector<Tensor<1,dim> > > &) const
{
  Assert (false, ExcPureFunctionCalled());
}


//////////////////////////////////////////////////////////////////////
// TensorFunction
//////////////////////////////////////////////////////////////////////

template <int rank_, int dim>
TensorFunction<rank_, dim>::TensorFunction (const double initial_time)
		:
		VectorFunction<dim>(pow(dim,rank_), initial_time)
{};



template <int rank_, int dim>
TensorFunction<rank_, dim>::~TensorFunction ()
{};



// template <int rank_, int dim>
// double
// TensorFunction<rank_, dim>::operator () (TensorIndex<rank_> i,
// 					 const Point<dim> &) const
// {
//   int k=i(0);
//   k++;
  
//   Assert (false, ExcPureFunctionCalled());
//   return 0;
// };


template <int rank_, int dim>
Tensor<rank_,dim>
TensorFunction<rank_, dim>::operator() (const Point<dim> &) const
{
  Assert (false, ExcPureFunctionCalled());
  return Tensor<rank_,dim>();
};


template <int rank_, int dim>
void
TensorFunction<rank_, dim>::value_list (const vector<Point<dim> > &points,
				 vector<Tensor<rank_,dim> > &values) const
{
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    values[i]  = this->operator() (points[i]);
};




template <int rank_, int dim>
Tensor<rank_+1,dim>
TensorFunction<rank_, dim>::gradient (const Point<dim> &) const
{
  Assert (false, ExcPureFunctionCalled());
  return Tensor<rank_+1,dim>();
};



template <int rank_, int dim>
void
TensorFunction<rank_, dim>::gradient_list (const vector<Point<dim> > &points,
				    vector<Tensor<rank_+1,dim> > &gradients) const
{
  Assert (gradients.size() == points.size(),
	  ExcVectorHasWrongSize(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i] = gradient(points[i]);
};

/*
template <int rank_, int dim> void
TensorFunction<rank_, dim>::value (const Point<dim>  &p, vector<double> &erg) const
{
  Tensor<rank_,dim> h = operator()(p);
  h.unroll(erg);
}
*/

template <int rank_, int dim> void
TensorFunction<rank_, dim>::value_list (const vector<Point<dim> > & points,
				 vector<vector<double> > & values) const
{
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    operator() (points[i]).unroll(values[i]);
  
}


template <int rank_, int dim> void
TensorFunction<rank_, dim>::gradient_list (const vector<Point<dim> > &,
				    vector<vector<Tensor<1,dim> > > &) const
{
  Assert (false, ExcPureFunctionCalled());
}



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
