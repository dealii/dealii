//----------------------------  block_vector.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_vector.templates.h  ---------------------------
#ifndef __deal2__block_vector_templates_h
#define __deal2__block_vector_templates_h


#include <base/config.h>
#include <base/memory_consumption.h>
#include <lac/block_vector.h>
#include <cmath>
#include <algorithm>


template <typename Number>
BlockVector<Number>::BlockVector (unsigned int n_blocks,
				  unsigned int block_size)
{
  reinit (n_blocks, block_size);
}



// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
  using namespace std;
#endif

template <typename Number>
BlockVector<Number>::BlockVector (const std::vector<unsigned int> &n)
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
void BlockVector<Number>::reinit (const unsigned int n_bl,
				  const unsigned int bl_sz,
				  const bool         fast)
{
  std::vector<unsigned int> n(n_bl, bl_sz);
  reinit(n, fast);
}


template <typename Number>
void BlockVector<Number>::reinit (const std::vector<unsigned int> &n,
				  const bool                       fast)
{
  block_indices.reinit (n);
  num_blocks = n.size();
  if (components.size() != num_blocks)
    components.resize(num_blocks);
  
  for (unsigned int i=0; i<num_blocks; ++i)
    components[i].reinit(n[i], fast);
}


template <typename Number>
template <typename Number2>
void BlockVector<Number>::reinit (const BlockVector<Number2>& v,
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
typename BlockVector<Number>::iterator
BlockVector<Number>::begin()
{
  return iterator(*this, 0U);
};



template <typename Number>
typename BlockVector<Number>::const_iterator
BlockVector<Number>::begin() const
{
  return const_iterator(*this, 0U);
};


template <typename Number>
typename BlockVector<Number>::iterator
BlockVector<Number>::end()
{
  return iterator(*this, size());
};



template <typename Number>
typename BlockVector<Number>::const_iterator
BlockVector<Number>::end() const
{
  return const_iterator(*this, size());
};



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
  return std::sqrt(norm_sqr());
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
template <typename Number2>
void BlockVector<Number>::equ (const Number a,
			       const BlockVector<Number2>& v)
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
BlockVector<Number>::operator = (const BlockVector<Number2>& v)
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
void BlockVector<Number>::print (std::ostream       &out,
				 const unsigned int  precision,
				 const bool          scientific,
				 const bool          across) const
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      if (across)
	out << 'C' << i << ':';
      else
	out << "Component " << i << std::endl;
      components[i].print(out, precision, scientific, across);
    }
}



template <typename Number>
void BlockVector<Number>::block_write (std::ostream &out) const
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].block_write(out);
    }
}



template <typename Number>
void BlockVector<Number>::block_read (std::istream &in)
{
  for (unsigned int i=0;i<num_blocks;++i)
    {
      components[i].block_read(in);
    }  
}


template <typename Number>
unsigned int
BlockVector<Number>::memory_consumption () const
{
  unsigned int mem = sizeof(num_blocks);
  for (unsigned int i=0; i<components.size(); ++i)
    mem += MemoryConsumption::memory_consumption (components[i]);
  mem += MemoryConsumption::memory_consumption (block_indices);
  return mem;
};



namespace BlockVectorIterators
{
  

  template <typename number, bool constness>
  Iterator<number,constness>::
  Iterator (BlockVectorType &parent,
	    const unsigned   global_index)
		  :
		  parent (&parent),
		  global_index (global_index)
  {
				     // find which block we are
				     // in. for this, take into
				     // account that it happens at
				     // times that people want to
				     // initialize iterators
				     // past-the-end
    if (global_index < parent.size())
      {
	const std::pair<unsigned int, unsigned int>
	  indices = parent.block_indices.global_to_local(global_index);
	current_block      = indices.first;
	index_within_block = indices.second;
	
	next_break_backward
	  = parent.block_indices.local_to_global (current_block, 0);
	next_break_forward
	  = (parent.block_indices.local_to_global (current_block, 0)
	     +parent.block_indices.block_size(current_block)-1);
      }
    else
				       // past the end. only have one
				       // value for this
      {
	this->global_index  = parent.size ();
	current_block       = parent.n_blocks();
	index_within_block  = 0;
	next_break_backward = global_index;
	next_break_forward  = static_cast<unsigned int>(-1);
      };
  };

  

  template <typename number, bool constness>
  void
  Iterator<number,constness>::move_forward ()
  {
    if (global_index != next_break_forward)
      ++index_within_block;
    else
      {
					 // ok, we traverse a boundary
					 // between blocks:
	index_within_block = 0;
	++current_block;

					 // break backwards is now old
					 // break forward
	next_break_backward = next_break_forward+1;

					 // compute new break forward
	if (current_block < parent->block_indices.size())
	  next_break_forward
	    += parent->block_indices.block_size(current_block);
	else
					   // if we are beyond the end,
					   // then move the next
					   // boundary arbitrarily far
					   // away
	  next_break_forward = static_cast<unsigned int>(-1);
      };
  
    ++global_index;
  };



  template <typename number, bool constness>
  void
  Iterator<number,constness>::move_backward ()
  {
    if (global_index != next_break_backward)
      --index_within_block;
    else
      if (current_block != 0)
	{
					   // ok, we traverse a boundary
					   // between blocks:
	  --current_block;
	  index_within_block = parent->block_indices.block_size(current_block)-1;
	
					   // break forwards is now old
					   // break backward
	  next_break_forward = next_break_backward-1;
	
					   // compute new break forward
	  next_break_backward
	    -= parent->block_indices.block_size (current_block);
	}
      else
					 // current block was 0, we now
					 // get into unspecified terrain
	{
	  --current_block;
	  index_within_block = static_cast<unsigned int>(-1);
	  next_break_forward = 0;
	  next_break_backward = 0;
	};
  
    --global_index;
  }; 
};


#endif
