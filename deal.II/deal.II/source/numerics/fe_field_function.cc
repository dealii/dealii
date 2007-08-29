//---------------------------------------------------------------------------
//    $Id: function_parser.h 14594 2007-03-22 20:17:41Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <numerics/fe_field_function.templates.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <hp/dof_handler.h>

#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{

  template class FEFieldFunction<deal_II_dimension, 
				 DoFHandler<deal_II_dimension>, 
				 Vector<double> >;

  template class FEFieldFunction<deal_II_dimension, 
				 DoFHandler<deal_II_dimension>, 
				 BlockVector<double> >;

  template class FEFieldFunction<deal_II_dimension, 
				 MGDoFHandler<deal_II_dimension>, 
				 Vector<double> >;

  template class FEFieldFunction<deal_II_dimension, 
				 MGDoFHandler<deal_II_dimension>, 
				 BlockVector<double> >;

#ifdef DEAL_II_USE_PETSC

  template class FEFieldFunction<deal_II_dimension, 
				 DoFHandler<deal_II_dimension>, 
				 PETScWrappers::Vector >;

  template class FEFieldFunction<deal_II_dimension, 
				 DoFHandler<deal_II_dimension>, 
				 PETScWrappers::BlockVector >;

  template class FEFieldFunction<deal_II_dimension, 
				 MGDoFHandler<deal_II_dimension>, 
				 PETScWrappers::Vector >;

  template class FEFieldFunction<deal_II_dimension, 
				 MGDoFHandler<deal_II_dimension>, 
				 PETScWrappers::BlockVector >;

#endif
  
}

DEAL_II_NAMESPACE_CLOSE
