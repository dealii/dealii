//----------------------------  log.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  log.cc  ---------------------------


#include <base/logstream.h>
#include <base/job_identifier.h>
#include <base/memory_consumption.h>

#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <iomanip>
#include <fstream>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif

LogStream deallog;


LogStream::LogStream()
		:
		std_out(&std::cerr), file(0), was_endl(true),
		std_depth(10000), file_depth(10000),
		print_utime(false), diff_utime(false),
		last_time (0.)
{
  prefixes.push("DEAL:");
  std_out->setf(std::ios::showpoint | std::ios::left);
}


void
LogStream::attach(std::ostream& o)
{
  file = &o;
  o.setf(std::ios::showpoint | std::ios::left);
  o << dealjobid();
}


void LogStream::detach ()
{
  file = 0;
};


std::ostream&
LogStream::get_console()
{
  return *std_out;
}


std::ostream&
LogStream::get_file_stream()
{
  Assert(file, ExcNoFileStreamGiven());
  return *file;
}


const std::string&
LogStream::get_prefix() const
{
  return prefixes.top();
}


void
LogStream::push (const std::string& text)
{
  std::string pre=prefixes.top();
  pre += text;
  pre += std::string(":");
  prefixes.push(pre);
};


void LogStream::pop ()
{
  if (prefixes.size() > 1)
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
LogStream::log_execution_time (bool flag)
{
  print_utime = flag;
}


void
LogStream::log_time_differences (bool flag)
{
  diff_utime = flag;
}


void
LogStream::print_line_head()
{
  rusage usage;
  double utime = 0.;
  if (print_utime)
    {
      getrusage(RUSAGE_SELF, &usage);
      utime = usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;
      if (diff_utime)
	{
	  double diff = utime - last_time;
	  last_time = utime;
	  utime = diff;
	}
    }

/*
 * The following lines were used for debugging a memory leak.
 * They work on Linux, not on Solaris, since the /proc filesystem
 * on Solaris is quite cryptic.
 *
 * Unfortunately, the information in /proc/pid/stat is updated slowly,
 * therefore, the information is quite unreliable.
 *
 * Furthermore, the cunstructor of ifstream caused another memory leak.
 *
 * Still, this code might be usefull sometimes, so I kept it here.
 * When we have more information about the kernel, this should be
 * incorporated properly. Suggestions are welcome!
 */
  
#ifdef DEALII_MEMORY_DEBUG
  static const pid_t id = getpid();
  
#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream statname;
#else
  std::ostrstream statname;
#endif
  
  statname << "/proc/" << id << "/stat" << std::ends;
  static long size;
  static string dummy;
  ifstream stat(statname.str());
				   // ignore 22 values
  stat >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >>
    dummy >> dummy >> dummy >> dummy >> dummy >>
    dummy >> dummy >> dummy >> dummy >> dummy >> dummy >>
    dummy >> dummy >> dummy >> dummy >> dummy >> size;
#endif
  
  const std::string& head = get_prefix();

  if (prefixes.size() <= std_depth)
    {
      if (print_utime)
	{
	  int p = std_out->width(5);
	  *std_out << utime << ':';
#ifdef DEALII_MEMORY_DEBUG
	  *std_out << size << ':';
#endif
	  std_out->width(p);
	}
      *std_out <<  head << ':';
    }
  
  if (file && (prefixes.size() <= file_depth))
    {
      if (print_utime)
	{
	  int p = file->width(6);
	  *file << utime << ':';
#ifdef DEALII_MEMORY_DEBUG
	  *file << size << ':';
#endif
	  file->width(p);
	}  
      *file << head << ':';
    }
}


unsigned int
LogStream::memory_consumption () const
{
  unsigned int mem = sizeof(*this);
				   // to determine size of stack
				   // elements, we have to copy the
				   // stack since we can't access
				   // elements from further below
  std::stack<std::string> tmp;
  while (tmp.size() > 0)
    {
      mem += MemoryConsumption::memory_consumption (tmp.top());
      tmp.pop ();
    };
  
  return mem;
};
