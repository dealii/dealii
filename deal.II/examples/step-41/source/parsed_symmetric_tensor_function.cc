/*
 * Immersed Boundary Problem:
 * 
 * Header Files
 *
 * Author: 
 * Luca Heltai <heltai@dimat.unipv.it>
 * =============
 * License: GPL.
 * =============
 * $Id: ibm_rhs.cc,v 1.34 2005/04/05 14:46:49 luca Exp $
 *
 */
#include "../include/parsed_symmetric_tensor_function.templates.h"

#if deal_II_dimension == 1
template <>
const SymmetricTensor<2, 1> & 
  ParsedSymmetricTensorFunction<2, 1>::operator()(const Point<1> &p) const
{ 
  t[0][0] = f.value(p,0);
  return t;
}


template <>
const SymmetricTensor<4, 1> & 
  ParsedSymmetricTensorFunction<4, 1>::operator()(const Point<1> &p) const
{ 
  t[0][0][0][0] = f.value(p, 0);
  return t;
}

#endif

#if deal_II_dimension == 2
template <>
const SymmetricTensor<2, 2> & 
  ParsedSymmetricTensorFunction<2, 2>::operator()(const Point<2> &p) const
{ 
  t[0][0] = f.value(p,0);
  t[0][1] = f.value(p,1);
  t[1][1] = f.value(p,2);
  return t;
}


template <>
const SymmetricTensor<4, 2> & 
  ParsedSymmetricTensorFunction<4, 2>::operator()(const Point<2> &p) const
{ 
  t[0][0][0][0] = f.value(p, 0);
  
  t[0][0][1][1] = f.value(p, 1);
  t[1][1][0][0] = f.value(p, 1);
  
  t[0][0][0][1] = f.value(p, 2);
  t[0][1][0][0] = f.value(p, 2);

  t[1][1][1][1] = f.value(p, 3);
  
  t[1][1][0][1] = f.value(p, 4);
  t[0][1][1][1] = f.value(p, 4);
  
  t[0][1][0][1] = f.value(p, 5);
  t[1][0][0][1] = f.value(p, 5);
  return t;
}

#endif

#if deal_II_dimension == 3

template <>
const SymmetricTensor<2, 3> & 
  ParsedSymmetricTensorFunction<2, 3>::operator()(const Point<3> &p) const
{ 
  t[0][0] = f.value(p,0);
  t[0][1] = f.value(p,1);
  t[0][2] = f.value(p,2);
  
  t[1][1] = f.value(p,3);
  t[1][2] = f.value(p,4);

  t[2][2] = f.value(p,5);

  return t;
}


template <>
const SymmetricTensor<4, 3> & 
  ParsedSymmetricTensorFunction<4, 3>::operator()(const Point<3> &p) const
{ 
  t[0][0][0][0] = f.value(p, 0);
  
  t[0][0][1][1] = f.value(p, 1);
  t[1][1][0][0] = f.value(p, 1);
  
  t[0][0][2][2] = f.value(p, 2);
  t[2][2][0][0] = f.value(p, 2);
  
  t[0][0][0][1] = f.value(p, 3);
  t[0][1][0][0] = f.value(p, 3);
  
  t[0][0][1][2] = f.value(p, 4);
  t[1][2][0][0] = f.value(p, 4);
  
  t[0][0][2][0] = f.value(p, 5);
  t[2][0][0][0] = f.value(p, 5);
  
  t[1][1][1][1] = f.value(p, 6);
  
  t[1][1][2][2] = f.value(p, 7);
  t[2][2][1][1] = f.value(p, 7);
  
  t[1][1][0][1] = f.value(p, 8);
  t[0][1][1][1] = f.value(p, 8);
  
  t[1][1][1][2] = f.value(p, 9);
  t[1][2][1][1] = f.value(p, 9);
  
  t[1][1][2][0] = f.value(p, 10);
  t[2][0][1][1] = f.value(p, 10);
  
  t[2][2][2][2] = f.value(p, 11);
  
  t[2][2][0][1] = f.value(p, 12);
  t[0][1][2][2] = f.value(p, 12);
  
  t[2][2][1][2] = f.value(p, 13);
  t[1][2][2][2] = f.value(p, 13);
  
  t[2][2][0][2] = f.value(p, 14);
  t[0][2][2][2] = f.value(p, 14);
  
  t[0][1][0][1] = f.value(p, 15);
  
  t[0][1][1][2] = f.value(p, 16);
  t[1][2][0][1] = f.value(p, 16);
  
  t[0][1][0][2] = f.value(p, 17);
  t[0][2][0][1] = f.value(p, 17);
  
  t[1][2][1][2] = f.value(p, 18);
  
  t[1][2][0][2] = f.value(p, 19);
  t[0][2][1][2] = f.value(p, 19);
  
  t[0][2][0][2] = f.value(p, 20);

  return t;
}

#endif

template class ParsedSymmetricTensorFunction<2,deal_II_dimension>;
template class ParsedSymmetricTensorFunction<4,deal_II_dimension>;
