// $Id$


#include <base/tensorfunction.h>
#include <vector>
#include <base/tensor.h>
#include <cmath>
#include <lac/vector.h>

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

template <int dim>
void
VectorFunction<dim>::value (const Point<dim>  &, Vector<double> &) const
{
  Assert (false, ExcPureFunctionCalled());
}


template <int dim>
void
VectorFunction<dim>::value_list (const vector<Point<dim> > &ps,
				 vector<Vector<double> > &us) const
{
  for (unsigned int i=0 ; i<ps.size() ; ++i)
    value(ps[i], us[i]);
}


template <int dim>
void
VectorFunction<dim>::gradient_list (const vector<Point<dim> > &,
				    vector<vector<Tensor<1,dim> > > &) const
{
  Assert (false, ExcPureFunctionCalled());
}

template <int dim>
VectorFunction<dim>::Extractor::Extractor(const VectorFunction<dim>& f,
					  unsigned int index)
		:
		vectorfunction(f),
		index(index)
{}

template <int dim>
double
VectorFunction<dim>::Extractor::operator() (const Point<dim>& p) const
{
  Vector<double> v(vectorfunction->n_components);
  vectorfunction->value(p,v);
  return v(index);
}


template <int dim>
Tensor<1,dim>
VectorFunction<dim>::Extractor::gradient (const Point<dim>&) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1,dim>();
}

template <int dim>
void
VectorFunction<dim>::Extractor::value_list (const vector<Point<dim> > &points,
					    vector<double> &values) const
{
  vector<Vector<double> > v(values.size(),
			    Vector<double>(vectorfunction->n_components));
  vectorfunction->value_list(p,v);
  for (unsigned int i=0 ; i<values.size() ; ++i)
    values[i] = v[i](index);
}


template <int dim>
void
VectorFunction<dim>::Extractor::gradient_list (const vector<Point<dim> > &points,
					       vector<Tensor<1,dim> > &gradients) const
{
  vector<vector<Tensor<1,dim> > > v(values.size(),
			    vector<Tensor<1,dim> >(vectorfunction->n_components));
  vectorfunction->value_list(p,v);
  for (unsigned int i=0 ; i<values.size() ; ++i)
    values[i] = v[i][index];
}

//////////////////////////////////////////////////////////////////////
// TensorFunction
//////////////////////////////////////////////////////////////////////

template <int rank, int dim>
TensorFunction<rank, dim>::TensorFunction (const double initial_time)
		:
		VectorFunction<dim>(pow(dim,rank), initial_time)
{};



template <int rank, int dim>
TensorFunction<rank, dim>::~TensorFunction ()
{};



// template <int rank, int dim>
// double
// TensorFunction<rank, dim>::operator () (TensorIndex<rank> i,
// 					 const Point<dim> &) const
// {
//   int k=i(0);
//   k++;
  
//   Assert (false, ExcPureFunctionCalled());
//   return 0;
// };


template <int rank, int dim>
Tensor<rank,dim>
TensorFunction<rank, dim>::operator() (const Point<dim> &) const
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
    values[i]  = this->operator() (points[i]);
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


template <int rank, int dim> void
TensorFunction<rank, dim>::value (const Point<dim>  &p,
				   Vector<double> &erg) const
{
  Tensor<rank,dim> h = operator()(p);
  h.unroll(erg);
}


template <int rank, int dim> void
TensorFunction<rank, dim>::value_list (const vector<Point<dim> > & points,
				 vector<Vector<double> > & values) const
{
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    operator() (points[i]).unroll(values[i]);
  
}


template <int rank, int dim> void
TensorFunction<rank, dim>::gradient_list (const vector<Point<dim> > &,
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
