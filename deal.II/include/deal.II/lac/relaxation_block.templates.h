//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__relaxation_block_templates_h
#define __deal2__relaxation_block_templates_h

#include <lac/relaxation_block.h>
#include <lac/full_matrix.h>

DEAL_II_NAMESPACE_OPEN

template <class MATRIX, typename inverse_type>
inline
RelaxationBlock<MATRIX,inverse_type>::AdditionalData::AdditionalData ()
		:
		relaxation(1.),
		invert_diagonal(true),
		same_diagonal(false),
		inversion(PreconditionBlockBase<inverse_type>::gauss_jordan)
{}


template <class MATRIX, typename inverse_type>
inline
RelaxationBlock<MATRIX,inverse_type>::AdditionalData::AdditionalData (
  const BlockList& bl,
  const double relaxation,
  const bool invert_diagonal,
  const bool same_diagonal)
		:
		block_list(&bl),
		relaxation(relaxation),
		invert_diagonal(invert_diagonal),
		same_diagonal(same_diagonal),
		inversion(PreconditionBlockBase<inverse_type>::gauss_jordan)
{}


template <class MATRIX, typename inverse_type>
inline
void
RelaxationBlock<MATRIX,inverse_type>::initialize (
  const MATRIX& M,
  const AdditionalData parameters)
{
  Assert (parameters.invert_diagonal, ExcNotImplemented());
  
  clear();
//  Assert (M.m() == M.n(), ExcNotQuadratic());
  A = &M;
  additional_data = parameters;
  
  this->reinit(additional_data.block_list->size(), 0, additional_data.same_diagonal,
	       additional_data.inversion);
  
  if (parameters.invert_diagonal)
    invert_diagblocks();
}


template <class MATRIX, typename inverse_type>
inline
void
RelaxationBlock<MATRIX,inverse_type>::clear ()
{
  A = 0;
  additional_data.block_list = 0;
  PreconditionBlockBase<inverse_type>::clear ();
}


template <class MATRIX, typename inverse_type>
inline
void
RelaxationBlock<MATRIX,inverse_type>::invert_diagblocks ()
{
  const MATRIX &M=*A;
  FullMatrix<inverse_type> M_cell;

  if (this->same_diagonal())
    {
      Assert(false, ExcNotImplemented());
    }
  else
    {
      for (unsigned int block=0;block<additional_data.block_list->size();++block)
	{
	  const unsigned int bs = additional_data.block_list->block_size(block);
	  M_cell.reinit(bs, bs);

					   // Copy rows for this block
					   // into the matrix for the
					   // diagonal block
	  BlockList::const_iterator row = additional_data.block_list->begin(block);
	  BlockList::const_iterator end = additional_data.block_list->end(block);
	  for (unsigned int row_cell=0; row_cell<bs; ++row_cell, ++row)
	    {
//TODO:[GK] Optimize here
	      for (typename MATRIX::const_iterator entry = M.begin(*row);
		   entry != M.end(*row); ++entry)
		{
		  const unsigned int column = entry->column();
		  const unsigned int col_cell = additional_data.block_list->local_index(block, column);
		  if (col_cell != numbers::invalid_unsigned_int)
		    M_cell(row_cell, col_cell) = entry->value();
		}
	    }
					   // Now M_cell contains the
					   // diagonal block. Now
					   // store it and its
					   // inverse, if so requested.
	  if (this->store_diagonals())
	    {
	      this->diagonal(block).reinit(bs, bs);
	      this->diagonal(block) = M_cell;
	    }
	  switch(this->inversion)
	    {
	      case PreconditionBlockBase<inverse_type>::gauss_jordan:
		    this->inverse(block).reinit(bs, bs);
		    this->inverse(block).invert(M_cell);
		    break;
	      case PreconditionBlockBase<inverse_type>::householder:
		    this->inverse_householder(block).initialize(M_cell);
		    break;
	      case PreconditionBlockBase<inverse_type>::svd:
		    this->inverse_svd(block).reinit(bs, bs);
		    this->inverse_svd(block) = M_cell;
		    this->inverse_svd(block).compute_inverse_svd(0.);
		    break;
	      default:
		    Assert(false, ExcNotImplemented());
	    }
	}
    }
  this->inverses_computed(true);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
inline
void
RelaxationBlock<MATRIX,inverse_type>::do_step (
  Vector<number2>       &dst,
  const Vector<number2> &prev,
  const Vector<number2> &src,
  const bool backward) const
{
  Assert (additional_data.invert_diagonal, ExcNotImplemented());
  
  const MATRIX &M=*this->A;
  Vector<number2> b_cell, x_cell;
  
  const unsigned int n_blocks = additional_data.block_list->size();
  for (unsigned int bi=0;bi<n_blocks;++bi)
    {
      const unsigned int block = backward ? (n_blocks - bi - 1) : bi;
      const unsigned int bs = additional_data.block_list->block_size(block);
      
      b_cell.reinit(bs);
      x_cell.reinit(bs);
				       // Collect off-diagonal parts
      BlockList::const_iterator row = additional_data.block_list->begin(block);
      BlockList::const_iterator end = additional_data.block_list->end(block);
      for (unsigned int row_cell=0; row_cell<bs; ++row_cell, ++row)
	{
	  b_cell(row_cell) = src(*row);
	  for (typename MATRIX::const_iterator entry = M.begin(*row);
	       entry != M.end(*row); ++entry)
	    b_cell(row_cell) -= entry->value() * prev(entry->column());
	}
				       // Apply inverse diagonal
      this->inverse_vmult(block, x_cell, b_cell);
				       // Store in result vector
      row=additional_data.block_list->begin(block);
      for (unsigned int row_cell=0; row_cell<bs; ++row_cell, ++row)
	dst(*row) = prev(*row) + additional_data.relaxation * x_cell(row_cell);
    }
}


//----------------------------------------------------------------------//

template <class MATRIX, typename inverse_type>
template <typename number2>
void RelaxationBlockSOR<MATRIX,inverse_type>::step (
  Vector<number2>       &dst,
  const Vector<number2> &src) const
{
  this->do_step(dst, dst, src, false);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void RelaxationBlockSOR<MATRIX,inverse_type>::Tstep (
  Vector<number2>       &dst,
  const Vector<number2> &src) const
{
  this->do_step(dst, dst, src, true);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void RelaxationBlockSSOR<MATRIX,inverse_type>::step (
  Vector<number2>       &dst,
  const Vector<number2> &src) const
{
  this->do_step(dst, dst, src, true);
  this->do_step(dst, dst, src, false);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void RelaxationBlockSSOR<MATRIX,inverse_type>::Tstep (
  Vector<number2>       &dst,
  const Vector<number2> &src) const
{
  this->do_step(dst, dst, src, false);
  this->do_step(dst, dst, src, true);
}



DEAL_II_NAMESPACE_CLOSE


#endif
