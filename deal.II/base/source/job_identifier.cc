//----------------------------  job_identifier.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  job_identifier.cc  ---------------------------


#include <base/job_identifier.h>
#include <ctime>

#if HAVE_GET_HOSTNAME
# include <unistd.h>
#endif

JobIdentifier dealjobid;


JobIdentifier::JobIdentifier()
{
  time_t t = time(0);
  id = std::string(program_id());

#if HAVE_GETHOSTNAME
  char name[100];
  gethostname(name,99);
  id += std::string(name) + std::string(" ");
#endif

  id += std::string(ctime(&t));
}


const std::string
JobIdentifier::operator ()() const
{
  return id;
}


