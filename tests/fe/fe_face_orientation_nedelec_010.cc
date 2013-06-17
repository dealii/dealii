//----------------------------------------------------------------------
//    $Id: fe_face_orientation_nedelec_010.cc 26378 2012-09-14 12:58:04Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------
//
// face_orientation = false, face_flip = true, face_rotation = false

#include "fe_face_orientation_nedelec.h"



int
main()
{
  initlog (__FILE__);
  deallog.threshold_double (1.e-10);
  run (false, true, false);
}
