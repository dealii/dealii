// $Id$

#include <base/tensor.h>
#include <cmath>

template <int dim> void
Tensor<1,dim>::unroll( vector<double>& result) const
{
  Assert(false,
	 ExcWrongVectorSize(1,1));
}


template <int rank_, int dim> void
Tensor<rank_, dim>::unroll( vector<double>& result) const
{
  Assert(result.size()==pow(dim,rank_),
	 ExcWrongVectorSize(pow(dim,rank_), result.size()));

  unsigned index = 0;
  unroll_recursion(result,index);
}


template <int rank_, int dim> void
Tensor<rank_, dim>::unroll_recursion( vector<double>& result, unsigned& index) const
{
  for (unsigned i=0; i<dim; ++i)
    {
      operator[](i).unroll_recursion(result, index);
    }
    
}


template<int dim> void
Tensor<1,dim>::unroll_recursion( vector<double>& result, unsigned& index) const
{
  for (unsigned i=0; i<dim; ++i)
    {
      cerr << "[" << index << ',' << operator[](i) << ']' << endl;
      result[index++] = operator[](i);
    }
  
}


template class Tensor<1, 1>;
template class Tensor<1, 2>;
template class Tensor<1, 3>;
template class Tensor<2, 1>;
template class Tensor<2, 2>;
template class Tensor<2, 3>;
template class Tensor<3, 1>;
template class Tensor<3, 2>;
template class Tensor<3, 3>;
template class Tensor<4, 1>;
template class Tensor<4, 2>;
template class Tensor<4, 3>;
