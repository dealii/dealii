// $Id$

#include <base/logstream.h>
#include <base/jobidentifier.h>



LogStream deallog;



LogStream::LogStream()
		: std_out(&cerr), file(0), was_endl(true),
		  std_depth(10000), file_depth(10000)
{
  prefixes.push("DEAL:");
}




void
LogStream::attach(ostream& o)
{
  file = &o;
  o << dealjobid();
}



void LogStream::detach ()
{
  file = 0;
};



void LogStream::pop ()
{
  prefixes.pop();
};



void LogStream::depth_console(unsigned n)
{
  std_depth = n;
};



void LogStream::depth_file(unsigned n)
{
  file_depth = n;
};



LogStream&
LogStream::operator << (void (f)(LogStream &))
{
  f(*this);
  return *this;
}
