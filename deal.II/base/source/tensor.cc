//----------------------------  tensor.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor.cc  ---------------------------


#include <base/tensor.h>
#include <cmath>
#include <lac/vector.h>


template <int dim>
void
Tensor<1,dim>::unroll( Vector<double>& result) const
{
  Assert(result.size()==dim,
	 ExcWrongVectorSize(dim, result.size()));

  unsigned index = 0;
  unroll_recursion(result,index);
}



template <int rank_, int dim>
void
Tensor<rank_, dim>::unroll( Vector<double>& result) const
{
  Assert(result.size()==pow(dim,rank_),
	 ExcWrongVectorSize(static_cast<unsigned int>(pow(dim,rank_)),
			    result.size()));

  unsigned index = 0;
  unroll_recursion(result,index);
}



template <int rank_, int dim>
void
Tensor<rank_, dim>::unroll_recursion( Vector<double>& result, unsigned& index) const
{
  for (unsigned i=0; i<dim; ++i)
    {
      operator[](i).unroll_recursion(result, index);
    }
    
}



template<int dim>
void
Tensor<1,dim>::unroll_recursion( Vector<double>& result, unsigned& index) const
{
  for (unsigned i=0; i<dim; ++i)
    result(index++) = operator[](i);  
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
