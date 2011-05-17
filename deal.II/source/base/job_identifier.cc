//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/base/job_identifier.h>
#include <ctime>

#if HAVE_GETHOSTNAME
# include <unistd.h>
#endif

DEAL_II_NAMESPACE_OPEN


JobIdentifier dealjobid;


JobIdentifier::JobIdentifier()
{
  time_t t = std::time(0);
  id = std::string("JobId ");

#if defined(HAVE_GETHOSTNAME) && !defined(DEAL_II_BROKEN_GETHOSTNAME) && !defined(DEBUG)
  char name[100];
  gethostname(name,99);
  id += std::string(name) + std::string(" ");
#else
  id += std::string("unknown ");
#endif

  id += std::string(std::ctime(&t));
}


const std::string
JobIdentifier::operator ()() const
{
  return id;
}


std::string
JobIdentifier::base_name(const char* filename)
{
  std::string name(filename);
  std::string::size_type pos = name.find(".");
  name.erase(pos, name.size());
  pos = name.rfind("/");
  if (pos < name.size())
    name.erase(0,pos);
  return name;
}



DEAL_II_NAMESPACE_CLOSE
