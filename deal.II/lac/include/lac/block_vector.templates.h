//----------------------------  block_vector.templates.h  ---------------------------
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
//----------------------------  block_vector.templates.h  ---------------------------
#ifndef __deal2__block_vector_templates_h
#define __deal2__block_vector_templates_h


#include <lac/block_vector.h>
#include <cmath>
#include <algorithm>


template <typename Number>
BlockVector<Number>::BlockVector (unsigned int num_blocks)
  : components(num_blocks),
    block_indices(num_blocks),
    num_blocks(num_blocks)
{}



template <typename Number>
BlockVector<Number>::BlockVector (const vector<unsigned int> &n)
  : block_indices(num_blocks)
{
  reinit (n, false);
}


template <typename Number>
BlockVector<Number>::BlockVector (const BlockVector<Number>& v)
  : components(v.num_blocks),
    block_indices(v.block_indices),
    num_blocks(v.num_blocks)    
{
  for (unsigned int i=0; i<num_blocks; ++i)
    components[i] = v.components[i];
}


// see the .h file for why this function was disabled
//
// template <typename Number>
// template <typename OtherNumber>
// BlockVector<Number>::Vector (const BlockVector< OtherNumber>& v) :
// 		dim(v.size()),
// 		maxdim(v.size()),
// 		val(0)
// {
//   if (dim)
//     {
//       val = new Number[maxdim];
//       Assert (val != 0, ExcOutOfMemory());
//       copy (v.begin(), v.end(), begin());
//     }
// }


template <typename Number>
void BlockVector<Number>::reinit (const vector<unsigned int> &n,
				  const bool                  fast)
{
  block_indices.reinit (n);
  num_blocks = n.size();
  if (components.size() != num_blocks)
    components.resize(num_blocks);
  
  for (unsigned int i=0; i<num_blocks; ++i)
    components[i].reinit(n[i], fast);
}


template <typename Number>
void BlockVector<Number>::reinit (const BlockVector<Number>& v,
					   const bool fast)
{
  block_indices = v.block_indices;
  num_blocks = v.num_blocks;
  if (components.size() != num_blocks)
    components.resize(num_blocks);
  
  for (unsigned int i=0;i<num_blocks;++i)
    components[i].reinit(v.components[i], fast);
}



template <typename Number>
BlockVector<Number>::~BlockVector ()
{}



template <typename Number>
void BlockVector<Number>::swap (BlockVector<Number> &v)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  for (unsigned int i=0; i<num_blocks; ++i)
    ::swap (components[i], v.components[i]);
  ::swap (block_indices, v.block_indices);
};



template <typename Number>
void BlockVector<Number>::clear ()
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].clear();
    }
}


template <typename Number>
bool BlockVector<Number>::all_zero () const
{
  bool result = true;
  for (unsigned int i=0;i<num_blocks;++i)
    {
      result = result && components[i].all_zero();
    }
  return result;
}


template <typename Number>
Number BlockVector<Number>::operator * (const BlockVector<Number>& v) const
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  Number sum = 0.;
  for (unsigned int i=0;i<num_blocks;++i)
    {
      sum += components[i]*v.components[i];
    }
  return sum;
}


template <typename Number>
Number BlockVector<Number>::norm_sqr () const
{
  Number sum = 0.;
  for (unsigned int i=0;i<num_blocks;++i)
    {
      sum += components[i].norm_sqr();
    }
  return sum;
};


template <typename Number>
Number BlockVector<Number>::mean_value () const
{
  Number sum = 0.;
  for (unsigned int i=0;i<num_blocks;++i)
    {
      sum += components[i].mean_value() * components[i].size();
    }
  return sum/size();
}


template <typename Number>
Number BlockVector<Number>::l1_norm () const
{
  Number sum = 0.;
  for (unsigned int i=0;i<num_blocks;++i)
    {
      sum += components[i].l1_norm();
    }
  return sum;
}


template <typename Number>
Number BlockVector<Number>::l2_norm () const
{
  return sqrt(norm_sqr());
}


template <typename Number>
Number BlockVector<Number>::linfty_norm () const
{
  Number sum = 0.;
  for (unsigned int i=0;i<num_blocks;++i)
    {
      Number newval = components[i].linfty_norm();
      if (sum<newval)
	sum = newval;
    }
  return sum;
}


template <typename Number>
BlockVector<Number>&
BlockVector<Number>::operator += (const BlockVector<Number>& v)
{
  add (v);
  return *this;
}


template <typename Number>
BlockVector<Number>&
BlockVector<Number>::operator -= (const BlockVector<Number>& v)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i] -= v.components[i];
    }
  return *this;
}


template <typename Number>
void BlockVector<Number>::add (const Number a)
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].add(a);
    }
}


template <typename Number>
void BlockVector<Number>::add (const BlockVector<Number>& v)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].add(v.components[i]);
    }
}


template <typename Number>
void BlockVector<Number>::add (const Number a,
			       const BlockVector<Number>& v)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].add(a, v.components[i]);
    }
}


template <typename Number>
void BlockVector<Number>::add (const Number a,
			       const BlockVector<Number>& v,
			       const Number b,
			       const BlockVector<Number>& w)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  Assert (num_blocks == w.num_blocks,
	  ExcDimensionMismatch(num_blocks, w.num_blocks));
  
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].add(a, v.components[i], b, w.components[i]);
    }
}


template <typename Number>
void BlockVector<Number>::sadd (const Number x,
				const BlockVector<Number>& v)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].sadd(x, v.components[i]);
    }
}


template <typename Number>
void BlockVector<Number>::sadd (const Number x, const Number a,
				const BlockVector<Number>& v)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].sadd(x, a, v.components[i]);
    }
}


template <typename Number>
void BlockVector<Number>::sadd (const Number x, const Number a,
				const BlockVector<Number>& v,
				const Number b,
				const BlockVector<Number>& w)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  Assert (num_blocks == w.num_blocks,
	  ExcDimensionMismatch(num_blocks, w.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].sadd(x, a, v.components[i], b, w.components[i]);
    }
}


template <typename Number>
void BlockVector<Number>::sadd (const Number x, const Number a,
				const BlockVector<Number>& v,
				const Number b,
				const BlockVector<Number>& w,
				const Number c,
				const BlockVector<Number>& y)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));  
  Assert (num_blocks == w.num_blocks,
	  ExcDimensionMismatch(num_blocks, w.num_blocks));
  Assert (num_blocks == y.num_blocks,
	  ExcDimensionMismatch(num_blocks, y.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].sadd(x, a, v.components[i],
			 b, w.components[i], c, y.components[i]);
    }
}


template <typename Number>
void BlockVector<Number>::scale (const Number factor)
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].scale(factor);
    }
}


template <typename Number>
void BlockVector<Number>::equ (const Number a,
			       const BlockVector<Number>& v,
			       const Number b,
			       const BlockVector<Number>& w)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  Assert (num_blocks == w.num_blocks,
	  ExcDimensionMismatch(num_blocks, w.num_blocks));  
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].equ( a, v.components[i], b, w.components[i]);
    }
}


template <typename Number>
void BlockVector<Number>::equ (const Number a,
			       const BlockVector<Number>& v)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].equ( a, v.components[i]);
    }
}


template <typename Number>
BlockVector<Number>& BlockVector<Number>::operator = (const Number s)
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i] = s;
    }
  return *this;
}


template <typename Number>
BlockVector<Number>&
BlockVector<Number>::operator = (const BlockVector<Number>& v)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i] = v.components[i];
    }
  return *this;
}


template <typename Number>
template<typename Number2>
BlockVector<Number>&
BlockVector<Number>::operator = (const BlockVector< Number2>& v)
{
  Assert (num_blocks == v.num_blocks,
	  ExcDimensionMismatch(num_blocks, v.num_blocks));
  
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i] = v.components[i];
    }
  return *this;
}


template <typename Number>
void BlockVector<Number>::print (ostream &out,
				 unsigned int precision,
				 bool scientific,
				 bool across) const
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      if (across)
	out << 'C' << i << ':';
      else
	out << "Component " << i << endl;
      components[i].print(out, precision, scientific, across);
    }
}


template <typename Number>
void BlockVector<Number>::block_write (ostream &out) const
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].block_write(out);
    }
}


template <typename Number>
void BlockVector<Number>::block_read (istream &in)
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].block_read(in);
    }  
}

#endif
