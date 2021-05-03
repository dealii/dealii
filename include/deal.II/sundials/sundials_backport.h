//-----------------------------------------------------------
//
//    Copyright (C) 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

/*
 * The functions in this file are based on an implementation distributed within
 * the SUNDIALS package, see the license here:
 * https://computing.llnl.gov/projects/sundials/license.
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David J. Gardner, Carol S. Woodward, and
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------*/

#ifndef dealii_sundials_backport_h
#define dealii_sundials_backport_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <sundials/sundials_nvector.h>

DEAL_II_NAMESPACE_OPEN
namespace SUNDIALS
{
  namespace internal
  {
    N_Vector
    N_VNewEmpty()
    {
      N_Vector     v   = new _generic_N_Vector;
      N_Vector_Ops ops = new _generic_N_Vector_Ops;

      /* initialize operations to nullptr */

      /* constructors, destructors, and utility operations */
      ops->nvgetvectorid     = nullptr;
      ops->nvclone           = nullptr;
      ops->nvcloneempty      = nullptr;
      ops->nvdestroy         = nullptr;
      ops->nvspace           = nullptr;
      ops->nvgetarraypointer = nullptr;
      ops->nvsetarraypointer = nullptr;
#  if DEAL_II_SUNDIALS_VERSION_GTE(5, 0, 0)
      ops->nvgetcommunicator = nullptr;
      ops->nvgetlength       = nullptr;
#  endif

      /* standard vector operations */
      ops->nvlinearsum    = nullptr;
      ops->nvconst        = nullptr;
      ops->nvprod         = nullptr;
      ops->nvdiv          = nullptr;
      ops->nvscale        = nullptr;
      ops->nvabs          = nullptr;
      ops->nvinv          = nullptr;
      ops->nvaddconst     = nullptr;
      ops->nvdotprod      = nullptr;
      ops->nvmaxnorm      = nullptr;
      ops->nvwrmsnormmask = nullptr;
      ops->nvwrmsnorm     = nullptr;
      ops->nvmin          = nullptr;
      ops->nvwl2norm      = nullptr;
      ops->nvl1norm       = nullptr;
      ops->nvcompare      = nullptr;
      ops->nvinvtest      = nullptr;
      ops->nvconstrmask   = nullptr;
      ops->nvminquotient  = nullptr;

#  if DEAL_II_SUNDIALS_VERSION_GTE(5, 0, 0)
      /* fused vector operations (optional) */
      ops->nvlinearcombination = nullptr;
      ops->nvscaleaddmulti     = nullptr;
      ops->nvdotprodmulti      = nullptr;

      /* vector array operations (optional) */
      ops->nvlinearsumvectorarray         = nullptr;
      ops->nvscalevectorarray             = nullptr;
      ops->nvconstvectorarray             = nullptr;
      ops->nvwrmsnormvectorarray          = nullptr;
      ops->nvwrmsnormmaskvectorarray      = nullptr;
      ops->nvscaleaddmultivectorarray     = nullptr;
      ops->nvlinearcombinationvectorarray = nullptr;

      /* local reduction operations (optional) */
      ops->nvdotprodlocal     = nullptr;
      ops->nvmaxnormlocal     = nullptr;
      ops->nvminlocal         = nullptr;
      ops->nvl1normlocal      = nullptr;
      ops->nvinvtestlocal     = nullptr;
      ops->nvconstrmasklocal  = nullptr;
      ops->nvminquotientlocal = nullptr;
      ops->nvwsqrsumlocal     = nullptr;
      ops->nvwsqrsummasklocal = nullptr;

#    if DEAL_II_SUNDIALS_VERSION_GTE(5, 4, 0)
      /* XBraid interface operations */
      ops->nvbufsize   = nullptr;
      ops->nvbufpack   = nullptr;
      ops->nvbufunpack = nullptr;
#    endif

#    if DEAL_II_SUNDIALS_VERSION_GTE(5, 3, 0)
      /* debugging functions (called when SUNDIALS_DEBUG_PRINTVEC is defined) */
      ops->nvprint     = nullptr;
      ops->nvprintfile = nullptr;
#    endif
#  endif

      /* attach ops and initialize content to nullptr */
      v->ops     = ops;
      v->content = nullptr;

      return v;
    }



    void
    N_VFreeEmpty(N_Vector v)
    {
      if (v == nullptr)
        return;

      /* free non-nullptr ops structure */
      if (v->ops)
        delete v->ops;
      v->ops = nullptr;

      /* free overall N_Vector object and return */
      delete v;
    }



    int
    N_VCopyOps(N_Vector w, N_Vector v)
    {
      /* Check that ops structures exist */
      if (w == nullptr || v == nullptr)
        return (-1);
      if (w->ops == nullptr || v->ops == nullptr)
        return (-1);

      /* Copy ops from w to v */

      /* constructors, destructors, and utility operations */
      v->ops->nvgetvectorid     = w->ops->nvgetvectorid;
      v->ops->nvclone           = w->ops->nvclone;
      v->ops->nvcloneempty      = w->ops->nvcloneempty;
      v->ops->nvdestroy         = w->ops->nvdestroy;
      v->ops->nvspace           = w->ops->nvspace;
      v->ops->nvgetarraypointer = w->ops->nvgetarraypointer;
      v->ops->nvsetarraypointer = w->ops->nvsetarraypointer;
#  if DEAL_II_SUNDIALS_VERSION_GTE(5, 0, 0)
      v->ops->nvgetcommunicator = w->ops->nvgetcommunicator;
      v->ops->nvgetlength       = w->ops->nvgetlength;
#  endif

      /* standard vector operations */
      v->ops->nvlinearsum    = w->ops->nvlinearsum;
      v->ops->nvconst        = w->ops->nvconst;
      v->ops->nvprod         = w->ops->nvprod;
      v->ops->nvdiv          = w->ops->nvdiv;
      v->ops->nvscale        = w->ops->nvscale;
      v->ops->nvabs          = w->ops->nvabs;
      v->ops->nvinv          = w->ops->nvinv;
      v->ops->nvaddconst     = w->ops->nvaddconst;
      v->ops->nvdotprod      = w->ops->nvdotprod;
      v->ops->nvmaxnorm      = w->ops->nvmaxnorm;
      v->ops->nvwrmsnormmask = w->ops->nvwrmsnormmask;
      v->ops->nvwrmsnorm     = w->ops->nvwrmsnorm;
      v->ops->nvmin          = w->ops->nvmin;
      v->ops->nvwl2norm      = w->ops->nvwl2norm;
      v->ops->nvl1norm       = w->ops->nvl1norm;
      v->ops->nvcompare      = w->ops->nvcompare;
      v->ops->nvinvtest      = w->ops->nvinvtest;
      v->ops->nvconstrmask   = w->ops->nvconstrmask;
      v->ops->nvminquotient  = w->ops->nvminquotient;

#  if DEAL_II_SUNDIALS_VERSION_GTE(5, 0, 0)
      /* fused vector operations */
      v->ops->nvlinearcombination = w->ops->nvlinearcombination;
      v->ops->nvscaleaddmulti     = w->ops->nvscaleaddmulti;
      v->ops->nvdotprodmulti      = w->ops->nvdotprodmulti;

      /* vector array operations */
      v->ops->nvlinearsumvectorarray     = w->ops->nvlinearsumvectorarray;
      v->ops->nvscalevectorarray         = w->ops->nvscalevectorarray;
      v->ops->nvconstvectorarray         = w->ops->nvconstvectorarray;
      v->ops->nvwrmsnormvectorarray      = w->ops->nvwrmsnormvectorarray;
      v->ops->nvwrmsnormmaskvectorarray  = w->ops->nvwrmsnormmaskvectorarray;
      v->ops->nvscaleaddmultivectorarray = w->ops->nvscaleaddmultivectorarray;
      v->ops->nvlinearcombinationvectorarray =
        w->ops->nvlinearcombinationvectorarray;

      /* local reduction operations */
      v->ops->nvdotprodlocal     = w->ops->nvdotprodlocal;
      v->ops->nvmaxnormlocal     = w->ops->nvmaxnormlocal;
      v->ops->nvminlocal         = w->ops->nvminlocal;
      v->ops->nvl1normlocal      = w->ops->nvl1normlocal;
      v->ops->nvinvtestlocal     = w->ops->nvinvtestlocal;
      v->ops->nvconstrmasklocal  = w->ops->nvconstrmasklocal;
      v->ops->nvminquotientlocal = w->ops->nvminquotientlocal;
      v->ops->nvwsqrsumlocal     = w->ops->nvwsqrsumlocal;
      v->ops->nvwsqrsummasklocal = w->ops->nvwsqrsummasklocal;

#    if DEAL_II_SUNDIALS_VERSION_GTE(5, 4, 0)
      /* XBraid interface operations */
      v->ops->nvbufsize   = w->ops->nvbufsize;
      v->ops->nvbufpack   = w->ops->nvbufpack;
      v->ops->nvbufunpack = w->ops->nvbufunpack;
#    endif

#    if DEAL_II_SUNDIALS_VERSION_GTE(5, 3, 0)
      /* debugging functions (called when SUNDIALS_DEBUG_PRINTVEC is defined) */
      v->ops->nvprint     = w->ops->nvprint;
      v->ops->nvprintfile = w->ops->nvprintfile;
#    endif
#  endif

      return (0);
    }
  } // namespace internal
} // namespace SUNDIALS
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS
#endif // dealii_sundials_sunlinsol_newempty_h
