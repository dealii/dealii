//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/trilinos_precondition_block.h>

#ifdef DEAL_II_USE_TRILINOS

#  include <lac/trilinos_sparse_matrix.h>
#  include <lac/trilinos_block_sparse_matrix.h>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  PreconditionBlockBase::PreconditionBlockBase()
  {}



  PreconditionBlockBase::~PreconditionBlockBase()
  {}

} /* end of namespace TrilinosWrappers */

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
