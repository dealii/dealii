//----------------------------  rt_projection_01.cc  ---------------------------
//    rt_projection_01.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rt_projection_01.cc  ---------------------------

/*
 * Project the function [1,1] onto a deformed grid and see whether the RT
 * elements can represent it exactly.
 */



char logname[] = "rt_projection_01/output";
#include "deformed_projection.h"


void test ()
{
  FE_RaviartThomas<2> fe (0);
  check (fe);
}
