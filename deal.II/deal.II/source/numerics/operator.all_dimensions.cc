//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/parameter_handler.h>
#include <base/logstream.h>
#include <lac/vector_memory.h>

#include <numerics/operator.templates.h>
#include <numerics/newton.templates.h>
#include <numerics/theta_timestepping.templates.h>

#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>


DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
#include "operator.all_dimensions.inst"
}

DEAL_II_NAMESPACE_CLOSE
