// $Id$

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
