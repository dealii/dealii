//----------------------------  vector.templates.h  ---------------------------
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
//----------------------------  vector.templates.h  ---------------------------
#ifndef __deal2__matrix_lib_templates_h
#define __deal2__matrix_lib_templates_h

#include <lac/matrix_lib.h>
#include <lac/vector.h>
#include <lac/block_vector.h>


template <typename number>
void
MeanValueFilter::vmult(Vector<number>& dst,
		       const Vector<number>& src) const
{
  Assert (dst.size() == src.size(),
	  ExcDimensionMismatch(dst.size(), src.size()));
  
  number mean = src.mean_value();

  for (unsigned int i=0;i<dst.size();++i)
    dst(i) = src(i) - mean;
}



template <typename number>
void
MeanValueFilter::vmult_add(Vector<number>& dst,
			   const Vector<number>& src) const
{
  Assert (dst.size() == src.size(),
	  ExcDimensionMismatch(dst.size(), src.size()));
  
  number mean = src.mean_value();
  
  for (unsigned int i=0;i<dst.size();++i)
    dst(i) += src(i) - mean;
}



template <typename number>
void
MeanValueFilter::vmult(BlockVector<number>& dst,
			   const BlockVector<number>& src) const
{
  Assert (component != static_cast<unsigned int>(-1),
	  ExcNotInitialized());
  
  Assert (dst.n_blocks() == src.n_blocks(),
	  ExcDimensionMismatch(dst.n_blocks(), src.n_blocks()));
  
  for (unsigned int i=0;i<dst.n_blocks();++i)
    if (i == component)
      vmult(dst.block(i), src.block(i));
    else
      dst.block(i) = src.block(i);
}



template <typename number>
void
MeanValueFilter::vmult_add(BlockVector<number>& dst,
			   const BlockVector<number>& src) const
{
  Assert (component != static_cast<unsigned int>(-1),
	  ExcNotInitialized());
  
  Assert (dst.n_blocks() == src.n_blocks(),
	  ExcDimensionMismatch(dst.n_blocks(), src.n_blocks()));
  
  for (unsigned int i=0;i<dst.n_blocks();++i)
    if (i == component)
      vmult_add(dst.block(i), src.block(i));
    else
      dst.block(i).add(src.block(i));
}

#endif
