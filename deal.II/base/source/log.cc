//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/logstream.h>
#include <base/job_identifier.h>
#include <base/memory_consumption.h>
#include <base/thread_management.h>

// include sys/resource.h for rusage(). Mac OS X needs sys/time.h then
// as well (strange), so include that, too.
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>


// on SunOS 4.x, getrusage is stated in the man pages and exists, but
// is not declared in resource.h. declare it ourselves
#ifdef NO_HAVE_GETRUSAGE
extern "C" { 
  int getrusage(int who, struct rusage* ru);
}
#endif

// When modifying the prefix list, we need to lock it just in case
// another thread tries to do the same.
DEAL_II_NAMESPACE_OPEN

namespace 
{
  Threads::ThreadMutex log_lock;
  Threads::ThreadMutex write_lock;
}


LogStream deallog;


LogStream::LogStream()
		:
		std_out(&std::cerr), file(0),
		std_depth(10000), file_depth(10000),
		print_utime(false), diff_utime(false),
		last_time (0.), double_threshold(0.), old_cerr(0)
{
  prefixes.push("DEAL:");
  std_out->setf(std::ios::showpoint | std::ios::left);
}


LogStream::~LogStream()
{
  if (old_cerr)
    std::cerr.rdbuf(old_cerr);

				   // on some systems, destroying the
				   // outstreams objects of deallog
				   // triggers some sort of memory
				   // corruption, in particular when
				   // we also link with Trilinos;
				   // since this happens at the very
				   // end of the program, we take the
				   // liberty to simply not do it by
				   // putting that object into a
				   // deliberate memory leak and
				   // instead destroying an empty
				   // object
#ifdef DEAL_II_USE_TRILINOS  
  if (this == &deallog)
    (new stream_map_type())->swap (outstreams);
#endif
}


LogStream &
LogStream::operator<< (std::ostream& (*p) (std::ostream&))
{
				   // do the work that is common to
				   // the operator<< functions
  print (p);

				   // next check whether this is the
				   // <tt>endl</tt> manipulator, and if so
				   // set a flag
  std::ostream & (* const p_endl) (std::ostream&) = &std::endl;
  if (p == p_endl)
    {
      Threads::ThreadMutex::ScopedLock lock(write_lock);
      print_line_head();
      std::ostringstream& stream = get_stream();
      if (prefixes.size() <= std_depth)
	*std_out << stream.str();
      
      if (file && (prefixes.size() <= file_depth))
	*file << stream.str() << std::flush;
      
				       // Start a new string
      stream.str("");
    }
  return *this;
}


std::ostringstream&
LogStream::get_stream()
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  unsigned int id = Threads::this_thread_id();

  std_cxx0x::shared_ptr<std::ostringstream>& sptr = outstreams[id];
  if (sptr == 0)
    {
      sptr = std_cxx0x::shared_ptr<std::ostringstream> (new std::ostringstream());
      sptr->setf(std::ios::showpoint | std::ios::left);
    } 
  return *sptr;
}



void
LogStream::attach(std::ostream& o)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  file = &o;
  o.setf(std::ios::showpoint | std::ios::left);
  o << dealjobid();
}


void LogStream::detach ()
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  file = 0;
}


void LogStream::log_cerr ()
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  if (old_cerr == 0)
    {
      old_cerr = std::cerr.rdbuf(file->rdbuf());
    } else {
      std::cerr.rdbuf(old_cerr);
      old_cerr = 0;
    }
}


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


bool
LogStream::has_file() const
{
  return (file != 0);
}


const std::string&
LogStream::get_prefix() const
{
  return prefixes.top();
}


void
LogStream::push (const std::string& text)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  std::string pre=prefixes.top();
  pre += text;
  pre += std::string(":");
  prefixes.push(pre);
}


void LogStream::pop ()
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  if (prefixes.size() > 1)
    prefixes.pop();
}


unsigned int
LogStream::depth_console (const unsigned n)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  const unsigned int h = std_depth;
  std_depth = n;
  return h;
}


unsigned int
LogStream::depth_file (const unsigned n)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  const unsigned int h = file_depth;
  file_depth = n;
  return h;
}


void
LogStream::threshold_double (const double t)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  double_threshold = t;
}


bool
LogStream::log_execution_time (const bool flag)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  const bool h = print_utime;
  print_utime = flag;
  return h;
}


bool
LogStream::log_time_differences (const bool flag)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  const bool h = diff_utime;
  diff_utime = flag;
  return h;
}


bool
LogStream::log_thread_id (const bool flag)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  const bool h = print_thread_id;
  print_thread_id = flag;
  return h;
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
 * on Solaris is quite cryptic. For other systems, we don't know.
 *
 * Unfortunately, the information in /proc/pid/stat is updated slowly,
 * therefore, the information is quite unreliable.
 *
 * Furthermore, the constructor of ifstream caused another memory leak.
 *
 * Still, this code might be useful sometimes, so I kept it here.
 * When we have more information about the kernel, this should be
 * incorporated properly. Suggestions are welcome!
 */
  
#ifdef DEALII_MEMORY_DEBUG
  static const pid_t id = getpid();
  
  std::ostringstream statname;
  statname << "/proc/" << id << "/stat";

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
  const unsigned int thread = Threads::this_thread_id();
  
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
      if (print_thread_id)
	*std_out << '[' << thread << ']';
      
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
      if (print_thread_id)
	*file << '[' << thread << ']';
      
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
}

DEAL_II_NAMESPACE_CLOSE
