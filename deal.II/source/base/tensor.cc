//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/tensor.h>
#include <cmath>
#include <lac/vector.h>

DEAL_II_NAMESPACE_OPEN


// storage for static variables
template <int dim, typename Number>
const unsigned int Tensor<1,dim,Number>::dimension;

template <int rank, int dim, typename Number>
const unsigned int Tensor<rank,dim,Number>::dimension;


template <int dim, typename Number>
void
Tensor<1,dim,Number>::unroll (Vector<Number> &result) const
{
  Assert (result.size()==dim,
          ExcDimensionMismatch(dim, result.size()));

  unsigned index = 0;
  unroll_recursion (result,index);
}



template <int rank_, int dim, typename Number>
void
Tensor<rank_, dim, Number>::unroll (Vector<Number> &result) const
{
  Assert(result.size()==std::pow(static_cast<Number>(dim),rank_),
	 ExcDimensionMismatch(static_cast<unsigned int>(std::pow(static_cast<Number>(dim),rank_)),
                              result.size()));

  unsigned index = 0;
  unroll_recursion (result, index);
}



template <int rank_, int dim, typename Number>
void
Tensor<rank_, dim, Number>::unroll_recursion (Vector<Number> &result,
					      unsigned int   &index) const
{
  for (unsigned i=0; i<dim; ++i)
    {
      operator[](i).unroll_recursion(result, index);
    }
}



template<int dim, typename Number>
void
Tensor<1,dim,Number>::unroll_recursion (Vector<Number> &result,
					unsigned int   &index) const
{
  for (unsigned i=0; i<dim; ++i)
    result(index++) = operator[](i);
}


template class Tensor<1, 1, double>;
template class Tensor<1, 2, double>;
template class Tensor<1, 3, double>;
template class Tensor<2, 1, double>;
template class Tensor<2, 2, double>;
template class Tensor<2, 3, double>;
template class Tensor<3, 1, double>;
template class Tensor<3, 2, double>;
template class Tensor<3, 3, double>;
template class Tensor<4, 1, double>;
template class Tensor<4, 2, double>;
template class Tensor<4, 3, double>;

DEAL_II_NAMESPACE_CLOSE
