//----------------------------  block_vector.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
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
BlockVector<Number>::BlockVector (const unsigned int n_blocks,
				  const unsigned int block_size)
{
  reinit (n_blocks, block_size);
}



template <typename Number>
BlockVector<Number>::BlockVector (const std::vector<unsigned int> &n)
{
  reinit (n, false);
}


template <typename Number>
BlockVector<Number>::BlockVector (const BlockVector<Number>& v)
                :
                BlockVectorBase<Vector<Number> > ()
{
  this->components.resize (v.num_blocks);
  this->block_indices = v.block_indices;
  num_blocks = v.num_blocks;
  
  for (unsigned int i=0; i<num_blocks; ++i)
    components[i] = v.components[i];
}


#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG    

template <typename Number>
template <typename OtherNumber>
BlockVector<Number>::BlockVector (const BlockVector<OtherNumber>& v)
{
  reinit (v, true);
  *this = v;
}

#endif

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
  this->block_indices.reinit (n);
  this->num_blocks = n.size();
  if (this->components.size() != num_blocks)
    this->components.resize(num_blocks);
  
  for (unsigned int i=0; i<this->num_blocks; ++i)
    this->components[i].reinit(n[i], fast);
}


template <typename Number>
template <typename Number2>
void BlockVector<Number>::reinit (const BlockVector<Number2>& v,
                                  const bool fast)
{
  this->block_indices = v.get_block_indices();
  this->num_blocks = v.n_blocks();
  if (this->components.size() != num_blocks)
    this->components.resize(num_blocks);
  
  for (unsigned int i=0;i<this->num_blocks;++i)
    this->block(i).reinit(v.block(i), fast);
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
}



#endif
