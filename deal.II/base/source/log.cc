//----------------------------  log.cc  ---------------------------
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
//----------------------------  log.cc  ---------------------------


#include <base/logstream.h>
#include <base/job_identifier.h>

#include <sys/resource.h>
#include <iomanip>

LogStream deallog;


LogStream::LogStream()
		: std_out(&cerr), file(0), was_endl(true),
		  std_depth(10000), file_depth(10000),
		  print_utime(false)
{
  prefixes.push("DEAL:");
  std_out->setf(ios::showpoint | ios::left);
}


void
LogStream::attach(ostream& o)
{
  file = &o;
  o.setf(ios::showpoint | ios::left);
  o << dealjobid();
}


void LogStream::detach ()
{
  file = 0;
};


ostream&
LogStream::get_console()
{
  return *std_out;
}


ostream&
LogStream::get_file_stream()
{
  Assert(file, ExcNoFileStreamGiven());
  return *file;
}


void
LogStream::push (const string& text)
{
  string pre=prefixes.top();
  pre += text;
  pre += string(":");
  prefixes.push(pre);
};


void LogStream::pop ()
{
  prefixes.pop();
};


void
LogStream::depth_console(unsigned n)
{
  std_depth = n;
};


void
LogStream::depth_file(unsigned n)
{
  file_depth = n;
};


void
LogStream::log_execution_time(bool flag)
{
  print_utime = flag;
}


void
LogStream::print_line_head()
{
  rusage usage;
  float utime = 0.;
  if (print_utime)
    {
      getrusage(RUSAGE_SELF, &usage);
      utime = usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;
    }
  
  if (prefixes.size() <= std_depth)
    {
      if (print_utime)
	{
	  int p = std_out->width(5);
	  *std_out << utime << ':';
	  std_out->width(p);
	}
    *std_out << prefixes.top() << ':';
    }
  
  if (file && (prefixes.size() <= file_depth))
    {
      if (print_utime)
	{
	  int p = file->width(6);
	  *file << utime << ':';
	  file->width(p);
	}  
    *file << prefixes.top() << ':';
    }
}


LogStream&
LogStream::operator << (void (f)(LogStream &))
{
  f(*this);
  return *this;
}


