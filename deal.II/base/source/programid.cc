//----------------------------  programid.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  programid.cc  ---------------------------


#include <base/job_identifier.h>

const char*
JobIdentifier::program_id()
{
  return "JobId ";
}
