// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/job_identifier.h>

#include <ctime>
#include <string>


#ifdef DEAL_II_HAVE_UNISTD_H
#  include <unistd.h>
#endif

DEAL_II_NAMESPACE_OPEN


const JobIdentifier &
JobIdentifier::get_dealjobid()
{
  static JobIdentifier dealjobid;
  return dealjobid;
}



JobIdentifier::JobIdentifier()
{
  std::time_t t = std::time(nullptr);
  id            = "JobId ";

#if defined(DEAL_II_HAVE_UNISTD_H) && defined(DEAL_II_HAVE_GETHOSTNAME)
  char name[100];
  gethostname(name, 99);
  id += std::string(name) + " ";
#else
  id += "unknown ";
#endif

  id += std::string(std::ctime(&t));
}


std::string
JobIdentifier::operator()() const
{
  return id;
}


std::string
JobIdentifier::base_name(const std::string &filename)
{
  std::string            name(filename);
  std::string::size_type pos;
  pos = name.rfind('/');
  if (pos != std::string::npos)
    name.erase(0, pos + 1);
  pos = name.rfind('.');
  if (pos != std::string::npos)
    name.erase(pos, name.size());
  return name;
}



DEAL_II_NAMESPACE_CLOSE
