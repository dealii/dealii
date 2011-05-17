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

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/algorithms/operator.templates.h>
#include <deal.II/algorithms/newton.templates.h>
#include <deal.II/algorithms/theta_timestepping.templates.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>


DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
#include "operator.inst"
}

DEAL_II_NAMESPACE_CLOSE
