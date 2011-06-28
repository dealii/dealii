//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__block_vector_templates_h
#define __deal2__block_vector_templates_h


#include <deal.II/base/config.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <cmath>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN

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
BlockVector<Number>::BlockVector (const BlockIndices& n)
{
  reinit (n, false);
}


template <typename Number>
BlockVector<Number>::BlockVector (const BlockVector<Number>& v)
                :
                BlockVectorBase<Vector<Number> > ()
{
  this->components.resize (v.n_blocks());
  this->block_indices = v.block_indices;
  
  for (unsigned int i=0; i<this->n_blocks(); ++i)
    this->components[i] = v.components[i];
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


#ifdef DEAL_II_USE_TRILINOS

template <typename Number>
BlockVector<Number>::BlockVector (const TrilinosWrappers::BlockVector &v)
{
  this->block_indices = v.get_block_indices();
  this->components.resize(this->n_blocks());

  for (unsigned int i=0; i<this->n_blocks(); ++i)
    this->components[i] = v.block(i);

  BaseClass::collect_sizes();
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
  if (this->components.size() != this->n_blocks())
    this->components.resize(this->n_blocks());
  
  for (unsigned int i=0; i<this->n_blocks(); ++i)
    this->components[i].reinit(n[i], fast);
}


template <typename Number>
void BlockVector<Number>::reinit (
  const BlockIndices& n,
  const bool fast)
{
  this->block_indices = n;
  if (this->components.size() != this->n_blocks())
    this->components.resize(this->n_blocks());
  
  for (unsigned int i=0; i<this->n_blocks(); ++i)
    this->components[i].reinit(n.block_size(i), fast);
}


template <typename Number>
template <typename Number2>
void BlockVector<Number>::reinit (const BlockVector<Number2>& v,
                                  const bool fast)
{
  this->block_indices = v.get_block_indices();
  if (this->components.size() != this->n_blocks())
    this->components.resize(this->n_blocks());
  
  for (unsigned int i=0;i<this->n_blocks();++i)
    this->block(i).reinit(v.block(i), fast);
}


template <typename Number>
BlockVector<Number>::~BlockVector ()
{}


#ifdef DEAL_II_USE_TRILINOS
template <typename Number>
inline
BlockVector<Number> &
BlockVector<Number>::operator = (const TrilinosWrappers::BlockVector &v)
{
  BaseClass::operator = (v);
  return *this;
}
#endif


template <typename Number>
void BlockVector<Number>::swap (BlockVector<Number> &v)
{
  Assert (this->n_blocks() == v.n_blocks(),
	  ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));
  
  for (unsigned int i=0; i<this->n_blocks(); ++i)
    dealii::swap (this->components[i], v.components[i]);
  dealii::swap (this->block_indices, v.block_indices);
}



template <typename Number>
void BlockVector<Number>::print (std::ostream       &out,
				 const unsigned int  precision,
				 const bool          scientific,
				 const bool          across) const
{
  for (unsigned int i=0;i<this->n_blocks();++i)
    {
      if (across)
	out << 'C' << i << ':';
      else
	out << "Component " << i << std::endl;
      this->components[i].print(out, precision, scientific, across);
    }
}



template <typename Number>
void BlockVector<Number>::block_write (std::ostream &out) const
{
  for (unsigned int i=0;i<this->n_blocks();++i)
    this->components[i].block_write(out);
}



template <typename Number>
void BlockVector<Number>::block_read (std::istream &in)
{
  for (unsigned int i=0;i<this->n_blocks();++i)
    this->components[i].block_read(in);
}



DEAL_II_NAMESPACE_CLOSE

#endif
