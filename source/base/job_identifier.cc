// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/job_identifier.h>
#include <ctime>

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#ifdef DEAM_II_MSVC
#  include <WinSock2.h>
#endif

DEAL_II_NAMESPACE_OPEN


JobIdentifier dealjobid;


JobIdentifier::JobIdentifier()
{
  time_t t = std::time(0);
  id = std::string("JobId ");

#if defined(HAVE_UNISTD_H) && defined(HAVE_GETHOSTNAME)
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
JobIdentifier::base_name(const char *filename)
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
