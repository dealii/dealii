//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__petsc_full_matrix_h
#define __deal2__petsc_full_matrix_h


#include <base/config.h>
#include <lac/exceptions.h>

#ifdef DEAL_II_USE_PETSC

#include <lac/petsc_matrix_base.h>


namespace PETScWrappers
{
/*! @addtogroup PETScWrappers
 *@{
 */

/**
 * Implementation of a sequential dense matrix class based on PETSC. All the
 * functionality is actually in the base class, except for the calls to
 * generate a sequential dense matrix. This is possible since PETSc only works
 * on an abstract matrix type and internally distributes to functions that do
 * the actual work depending on the actual matrix type (much like using
 * virtual functions). Only the functions creating a matrix of specific type
 * differ, and are implemented in this particular class.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class FullMatrix : public MatrixBase
  {
    public:
                                       /**
                                        * Create a full matrix of dimensions
                                        * @p m times @p n.
                                        */
      FullMatrix (const unsigned int m,
                  const unsigned int n);
  };
  
/*@}*/
}

#endif // DEAL_II_USE_PETSC

/*----------------------------   petsc_full_matrix.h     ---------------------------*/

#endif
/*----------------------------   petsc_full_matrix.h     ---------------------------*/
