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
  Threads::Mutex log_lock;
  Threads::Mutex write_lock;
}


// The standard log object of deal.II:
LogStream deallog;


LogStream::LogStream()
  :
  std_out(&std::cerr),
  file(0),
  std_depth(10000),
  file_depth(10000),
  print_utime(false),
  diff_utime(false),
  last_time (0.),
  double_threshold(0.),
  float_threshold(0.),
  offset(0),
  old_cerr(0),
  at_newline(true)
{
  get_prefixes().push("DEAL:");

#if defined(HAVE_UNISTD_H) && defined(HAVE_TIMES)
  reference_time_val = 1./sysconf(_SC_CLK_TCK) * times(&reference_tms);
#endif

}


LogStream::~LogStream()
{
  // if there was anything left in the stream that is current to this
  // thread, make sure we flush it before it gets lost
  {
    if (get_stream().str().length() > 0)
      {
        // except the situation is not quite that simple. if this object is
        // the 'deallog' object, then it is destroyed upon exit of the
        // program. since it's defined in a shared library that depends on
        // libstdc++.so, destruction happens before destruction of
        // std::cout/cerr, but after all file variables defined in user
        // programs have been destroyed. in other words, if we get here and
        // the object being destroyed is 'deallog' and if 'deallog' is
        // associated with a file stream, then we're in trouble: we'll try
        // to write to a file that doesn't exist any more, and we're likely
        // going to crash (this is tested by base/log_crash_01). rather
        // than letting it come to this, print a message to the screen
        // (note that we can't issue an assertion here either since Assert
        // may want to write to 'deallog' itself, and AssertThrow will
        // throw an exception that can't be caught)
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
}


void
LogStream::test_mode(bool on)
{
  Threads::Mutex::ScopedLock lock(log_lock);
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
LogStream::operator<< (std::ostream& (*p) (std::ostream &))
{

  std::ostringstream &stream = get_stream();

  // Print to the internal stringstream:
  stream << p;


  // This is a bloody hack until LogStream got reimplemented as a proper
  // child of std::streambuf (or similar).
  //
  // The problem is that at this point we would like to know whether an
  // std::flush or std::endl has called us, however, there is no way to
  // detect this in a sane manner.
  //
  // The obvious idea to compare function pointers,
  //   std::ostream & (* const p_flush) (std::ostream &) = &std::flush;
  //   p == p_flush ? ...,
  // is wrong as there doesn't has to be a _single_ std::flush instance...
  // there could be multiple of it. And in fact, LLVM's libc++ implements
  // std::flush and std::endl in a way that every shared library and
  // executable has its local copy... fun...
  //
  // - Maier, 2013

  class QueryStreambuf : public std::streambuf
  {
    // Implement a minimalistic stream buffer that only stores the fact
    // whether overflow or sync was called
  public:
    QueryStreambuf()
      : flushed_(false), newline_written_(false)
    {
    }
    bool flushed()
    {
      return flushed_;
    }
    bool newline_written()
    {
      return newline_written_;
    }
  private:
    int_type overflow(int_type ch)
    {
      newline_written_ = true;
      return ch;
    }
    int sync()
    {
      flushed_ = true;
      return 0;
    }
    bool flushed_;
    bool newline_written_;
  } query_streambuf;

  {
    // and initialize an ostream with this streambuf:
    std::ostream inject (&query_streambuf);
    inject << p;
  }

  if (query_streambuf.flushed())
    {
      Threads::Mutex::ScopedLock lock(write_lock);

      // Print the line head in case of a previous newline:
      if (at_newline)
        print_line_head();

      at_newline = query_streambuf.newline_written();

      if (get_prefixes().size() <= std_depth)
        *std_out << stream.str();

      if (file && (get_prefixes().size() <= file_depth))
        *file << stream.str() << std::flush;

      // Start a new string:
      stream.str("");
    }

  return *this;
}


void
LogStream::attach(std::ostream &o,
                  const bool    print_job_id)
{
  Threads::Mutex::ScopedLock lock(log_lock);
  file = &o;
  o.setf(std::ios::showpoint | std::ios::left);
  if (print_job_id)
    o << dealjobid();
}


void LogStream::detach ()
{
  Threads::Mutex::ScopedLock lock(log_lock);
  file = 0;
}


void LogStream::log_cerr ()
{
  Threads::Mutex::ScopedLock lock(log_lock);
  if (old_cerr == 0)
    {
      old_cerr = std::cerr.rdbuf(file->rdbuf());
    }
  else
    {
      std::cerr.rdbuf(old_cerr);
      old_cerr = 0;
    }
}


std::ostream &
LogStream::get_console()
{
  return *std_out;
}


std::ostream &
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


const std::string &
LogStream::get_prefix() const
{
  static std::string empty_string;

  if (get_prefixes().size() > 0)
    return get_prefixes().top();
  else
    return empty_string;
}


void
LogStream::push (const std::string &text)
{
  std::string pre;
  if (get_prefixes().size() > 0)
    pre = get_prefixes().top();

  pre += text;
  pre += std::string(":");
  get_prefixes().push(pre);
}


void LogStream::pop ()
{
  if (get_prefixes().size() > 0)
    get_prefixes().pop();
}


std::ios::fmtflags
LogStream::flags(const std::ios::fmtflags f)
{
  return get_stream().flags (f);
}


std::streamsize
LogStream::precision (const std::streamsize prec)
{
  return get_stream().precision (prec);
}


std::streamsize
LogStream::width (const std::streamsize wide)
{
  return get_stream().width (wide);
}


unsigned int
LogStream::depth_console (const unsigned int n)
{
  Threads::Mutex::ScopedLock lock(log_lock);
  const unsigned int h = std_depth;
  std_depth = n;
  return h;
}


unsigned int
LogStream::depth_file (const unsigned int n)
{
  Threads::Mutex::ScopedLock lock(log_lock);
  const unsigned int h = file_depth;
  file_depth = n;
  return h;
}


void
LogStream::threshold_double (const double t)
{
  Threads::Mutex::ScopedLock lock(log_lock);
  double_threshold = t;
}


void
LogStream::threshold_float (const float t)
{
  Threads::Mutex::ScopedLock lock(log_lock);
  float_threshold = t;
}


bool
LogStream::log_execution_time (const bool flag)
{
  Threads::Mutex::ScopedLock lock(log_lock);
  const bool h = print_utime;
  print_utime = flag;
  return h;
}


bool
LogStream::log_time_differences (const bool flag)
{
  Threads::Mutex::ScopedLock lock(log_lock);
  const bool h = diff_utime;
  diff_utime = flag;
  return h;
}


bool
LogStream::log_thread_id (const bool flag)
{
  Threads::Mutex::ScopedLock lock(log_lock);
  const bool h = print_thread_id;
  print_thread_id = flag;
  return h;
}

std::stack<std::string> &
LogStream::get_prefixes() const
{
#ifdef DEAL_II_WITH_THREADS
  bool exists = false;
  std::stack<std::string> &local_prefixes = prefixes.get(exists);

  // If this is a new locally stored stack, copy the "blessed" prefixes
  // from the initial thread that created logstream.
  if (! exists)
    {
      const tbb::enumerable_thread_specific<std::stack<std::string> > &impl
        = prefixes.get_implementation();

      // The thread that created this LogStream object should be the first
      // in tbb's enumerable_thread_specific containter.
      const tbb::enumerable_thread_specific<std::stack<std::string> >::const_iterator first_elem
        = impl.begin();

      if (first_elem != impl.end())
        {
          local_prefixes = *first_elem;
        }
    }

  return local_prefixes;

#else
  return prefixes.get();
#endif
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

  const std::string &head = get_prefix();
  const unsigned int thread = Threads::this_thread_id();

  if (get_prefixes().size() <= std_depth)
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

      if (head.size() > 0)
        *std_out <<  head << ':';
    }

  if (file && (get_prefixes().size() <= file_depth))
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

      if (head.size() > 0)
        *file << head << ':';
    }
}


void
LogStream::timestamp ()
{
  struct tms current_tms;
#if defined(HAVE_UNISTD_H) && defined(HAVE_TIMES)
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
  // TODO
  Assert(false, ExcNotImplemented());

  std::size_t mem = sizeof(*this);
  // to determine size of stack
  // elements, we have to copy the
  // stack since we can't access
  // elements from further below
//   std::stack<std::string> tmp = prefixes;
//   while (tmp.empty() == false)
//     {
//       mem += MemoryConsumption::memory_consumption (tmp.top());
//       tmp.pop ();
//     }

  return mem;
}

DEAL_II_NAMESPACE_CLOSE
