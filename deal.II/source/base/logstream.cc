//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2009, 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/base/logstream.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/thread_management.h>

#ifdef HAVE_SYS_RESOURCE_H
#  include <sys/resource.h>
#endif

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>


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
                last_time (0.), double_threshold(0.), float_threshold(0.),
                offset(0), old_cerr(0)
{
  prefixes.push("DEAL:");
  std_out->setf(std::ios::showpoint | std::ios::left);
#if defined(HAVE_TIMES) && defined(HAVE_UNISTD_H)
  reference_time_val = 1./sysconf(_SC_CLK_TCK) * times(&reference_tms);
#endif
}


LogStream::~LogStream()
{
                                   // if there was anything left in
                                   // the stream that is current to
                                   // this thread, make sure we flush
                                   // it before it gets lost
  {
    const unsigned int id = Threads::this_thread_id();
    if ((outstreams.find(id) != outstreams.end())
        &&
        (*outstreams[id] != 0)
        &&
        (outstreams[id]->str().length() > 0))
      {
					 // except the situation is
					 // not quite that simple. if
					 // this object is the
					 // 'deallog' object, then it
					 // is destroyed upon exit of
					 // the program. since it's
					 // defined in a shared
					 // library that depends on
					 // libstdc++.so, destruction
					 // happens before destruction
					 // of std::cout/cerr, but
					 // after all file variables
					 // defined in user programs
					 // have been destroyed. in
					 // other words, if we get
					 // here and the object being
					 // destroyed is 'deallog' and
					 // if 'deallog' is associated
					 // with a file stream, then
					 // we're in trouble: we'll
					 // try to write to a file
					 // that doesn't exist any
					 // more, and we're likely
					 // going to crash (this is
					 // tested by
					 // base/log_crash_01). rather
					 // than letting it come to
					 // this, print a message to
					 // the screen (note that we
					 // can't issue an assertion
					 // here either since Assert
					 // may want to write to
					 // 'deallog' itself, and
					 // AssertThrow will throw an
					 // exception that can't be
					 // caught)
	if ((this == &deallog) && file)
	  std::cerr << ("You still have content that was written to 'deallog' "
			"but not flushed to the screen or a file while the "
			"program is being terminated. This would lead to a "
			"segmentation fault. Make sure you flush the "
			"content of the 'deallog' object using 'std::endl' "
			"before the end of the program.")
		    << std::endl;
	else
	  *this << std::endl;
      }
  }

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
    {
      stream_map_type * dummy = new stream_map_type();
      dummy->swap (outstreams);
      delete dummy;
    }
#endif
}


void
LogStream::test_mode(bool on)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  if (on)
    {
      double_threshold = 1.e-10;
      float_threshold = 1.e-7;
      offset = 1.e-7;
    }
  else
    {
      double_threshold = 0.;
      float_threshold = 0.;
      offset = 0.;
    }
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
//TODO: use a ThreadLocalStorage object here
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  const unsigned int id = Threads::this_thread_id();

                                   // if necessary allocate a stream object
  if (outstreams.find (id) == outstreams.end())
    {
      outstreams[id].reset (new std::ostringstream());
      outstreams[id]->setf(std::ios::showpoint | std::ios::left);
    }
  return *outstreams[id];
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
LogStream::depth_console (const unsigned int n)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  const unsigned int h = std_depth;
  std_depth = n;
  return h;
}


unsigned int
LogStream::depth_file (const unsigned int n)
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


void
LogStream::threshold_float (const float t)
{
  Threads::ThreadMutex::ScopedLock lock(log_lock);
  float_threshold = t;
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
#ifdef HAVE_SYS_RESOURCE_H
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
#else
//TODO[BG]: Do something useful here
  double utime = 0.;
#endif

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


void
LogStream::timestamp ()
{
  struct tms current_tms;
#if defined(HAVE_TIMES) && defined(HAVE_UNISTD_H)
  const clock_t tick = sysconf(_SC_CLK_TCK);
  const double time = 1./tick * times(&current_tms);
#else
  const double time = 0.;
  const unsigned int tick = 100;
#endif
  (*this) << "Wall: " << time - reference_time_val
          << " User: " << 1./tick * (current_tms.tms_utime - reference_tms.tms_utime)
          << " System: " << 1./tick * (current_tms.tms_stime - reference_tms.tms_stime)
          << " Child-User: " << 1./tick * (current_tms.tms_cutime - reference_tms.tms_cutime)
          << " Child-System: " << 1./tick * (current_tms.tms_cstime - reference_tms.tms_cstime)
          << std::endl;
}


std::size_t
LogStream::memory_consumption () const
{
  std::size_t mem = sizeof(*this);
                                   // to determine size of stack
                                   // elements, we have to copy the
                                   // stack since we can't access
                                   // elements from further below
  std::stack<std::string> tmp = prefixes;
  while (tmp.empty() == false)
    {
      mem += MemoryConsumption::memory_consumption (tmp.top());
      tmp.pop ();
    }

  return mem;
}

DEAL_II_NAMESPACE_CLOSE
