//----------------------------  multigrid.all_dimensions.cc  ---------------------------
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
//----------------------------  multigrid.all_dimensions.cc  ---------------------------



#include <grid/tria.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <multigrid/multigrid.h>
#include <multigrid/multigrid.templates.h>
#include <multigrid/mg_smoother.h>
#include <lac/vector.h>


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


