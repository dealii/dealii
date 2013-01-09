//----------------------------  abfprojection_01.cc  ---------------------------
//    abfprojection_01.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  abfprojection_01.cc  ---------------------------

/*
 * Project the function [1,1] onto a deformed grid and see whether the
 * FESystem elements can represent it exactly. This shouldn't be a surprise,
 * but it is nice to compare with the RT and ABF elements
 */



char logname[] = "q_dg0_projection_01/output";
#include "deformed_projection.h"


void test ()
{
  FESystem<2> fe (FE_Q_DG0<2>(3), 2);
  check (fe);
}
