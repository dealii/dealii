//----------------------------  multigrid.all_dimensions.cc  ---------------------------
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
//----------------------------  multigrid.all_dimensions.cc  ---------------------------



#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <multigrid/mg_transfer.h>


MGTransferPrebuilt::~MGTransferPrebuilt () 
{};


void MGTransferPrebuilt::prolongate (const unsigned int   to_level,
				     Vector<double>       &dst,
				     const Vector<double> &src) const 
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[to_level-1].vmult (dst, src);
};


void MGTransferPrebuilt::restrict_and_add (const unsigned int   from_level,
					   Vector<double>       &dst,
					   const Vector<double> &src) const 
{
  Assert ((from_level >= 1) && (from_level<=prolongation_matrices.size()),
	  ExcIndexRange (from_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[from_level-1].Tvmult_add (dst, src);
};





MGTransferBlock::~MGTransferBlock () 
{};


void MGTransferBlock::prolongate (const unsigned int   to_level,
				     BlockVector<double>       &dst,
				     const BlockVector<double> &src) const 
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

  unsigned int k=0;
  for (unsigned int b=0; b<src.n_blocks();++b)
    {
      if (!selected[k])
	++k;
      prolongation_matrices[to_level-1].block(k,k).vmult (dst.block(b), src.block(b));
      ++k;
    }
};


void MGTransferBlock::restrict_and_add (const unsigned int   from_level,
					   BlockVector<double>       &dst,
					   const BlockVector<double> &src) const 
{
  Assert ((from_level >= 1) && (from_level<=prolongation_matrices.size()),
	  ExcIndexRange (from_level, 1, prolongation_matrices.size()+1));

  unsigned int k=0;
  for (unsigned int b=0; b<src.n_blocks();++b)
    {
      if (!selected[k])
	++k;
      prolongation_matrices[from_level-1].block(k,k).Tvmult_add (dst.block(b), src.block(b));
      ++k;
    }
};




MGTransferSelect::~MGTransferSelect () 
{};


void MGTransferSelect::prolongate (const unsigned int   to_level,
				     Vector<double>       &dst,
				     const Vector<double> &src) const 
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

      prolongation_matrices[to_level-1].block(selected, selected)
	.vmult (dst, src);
};


void MGTransferSelect::restrict_and_add (const unsigned int   from_level,
					   Vector<double>       &dst,
					   const Vector<double> &src) const 
{
  Assert ((from_level >= 1) && (from_level<=prolongation_matrices.size()),
	  ExcIndexRange (from_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[from_level-1].block(selected, selected)
    .Tvmult_add (dst, src);
};
