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

#include <lac/vector.h>
#include <lac/block_vector.h>

#ifdef DEAL_II_USE_PETSC
#  include <lac/petsc_vector.h>
#  include <lac/petsc_block_vector.h>
#endif

#ifdef DEAL_II_USE_TRILINOS
#  include <lac/trilinos_vector.h>
#  include <lac/trilinos_block_vector.h>
#endif

DEAL_II_NAMESPACE_OPEN

using namespace Algorithms;

#include "operators.inst"

DEAL_II_NAMESPACE_CLOSE
