// $Id$

#include <base/logstream.h>
#include <base/jobidentifier.h>

LogStream deallog;

LogStream::LogStream()
		: std(cerr), file(0), was_endl(true),
		  std_depth(10000), file_depth(10000)
{
  prefixes.push("DEAL:");
}

void
LogStream::attach(ostream& o)
{
  file = &o;
  o << "DEAL:" << dealjobid();
}
