//----------------------------  job_identifier.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  job_identifier.cc  ---------------------------


#include <base/job_identifier.h>
#include <ctime>


JobIdentifier dealjobid;

JobIdentifier::JobIdentifier()
{
  time_t t = time(0);
  id = string(program_id()) + string(ctime(&t));
}

JobIdentifier::~JobIdentifier()
{}

const string
JobIdentifier::operator ()() const
{
  return id;
}
