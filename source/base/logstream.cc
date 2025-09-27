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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/thread_management.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stack>
#include <thread>


DEAL_II_NAMESPACE_OPEN

namespace
{
  std::mutex log_lock;
  std::mutex write_lock;
} // namespace


// The standard log object of deal.II:
LogStream deallog;



LogStream::Prefix::Prefix(const std::string &text)
  : stream(&deallog)
{
  stream->push(text);
}



LogStream::Prefix::Prefix(const std::string &text, LogStream &s)
  : stream(&s)
{
  stream->push(text);
}



LogStream::Prefix::~Prefix()
{
  // destructors may not throw exceptions. if the pop() function
  // actually does throw one, then ignore it
  try
    {
      stream->pop();
    }
  catch (...)
    {
      AssertNothrow(false,
                    ExcMessage(
                      "An exception occurred in LogStream::Prefix::~Prefix."));
    }
}



LogStream::LogStream()
  : parent_thread(std::this_thread::get_id())
  , std_out(&std::cout)
  , file(nullptr)
  , std_depth(0)
  , file_depth(10000)
  , print_thread_id(false)
  , at_newline(true)
{
  get_prefixes().emplace("DEAL:");
}



LogStream::~LogStream()
{
  // if there was anything left in the stream that is current to this
  // thread, make sure we flush it before it gets lost
  {
    if (get_stream().str().size() > 0)
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
        if ((this == &deallog) && (file != nullptr))
          *std_out << ("You still have content that was written to 'deallog' "
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
}


LogStream &
LogStream::operator<<(std::ostream &(*p)(std::ostream &))
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
      : flushed_(false)
      , newline_written_(false)
    {}
    bool
    flushed()
    {
      return flushed_;
    }
    bool
    newline_written()
    {
      return newline_written_;
    }

  private:
    int_type
    overflow(int_type ch) override
    {
      newline_written_ = true;
      return ch;
    }
    int
    sync() override
    {
      flushed_ = true;
      return 0;
    }
    bool flushed_;
    bool newline_written_;
  } query_streambuf;

  {
    // and initialize an ostream with this streambuf:
    std::ostream inject(&query_streambuf);
    inject << p;
  }

  if (query_streambuf.flushed())
    {
      std::lock_guard<std::mutex> lock(write_lock);

      // Print the line head in case of a previous newline:
      if (at_newline)
        print_line_head();

      at_newline = query_streambuf.newline_written();

      if (get_prefixes().size() <= std_depth)
        *std_out << stream.str();

      if ((file != nullptr) && (get_prefixes().size() <= file_depth))
        *file << stream.str() << std::flush;

      // Start a new string:
      stream.str("");
    }

  return *this;
}


void
LogStream::attach(std::ostream                 &o,
                  const bool                    print_job_id,
                  const std::ios_base::fmtflags flags)
{
  std::lock_guard<std::mutex> lock(log_lock);
  file = &o;
  o.setf(flags);
  if (print_job_id)
    o << JobIdentifier::get_dealjobid()();
}


void
LogStream::detach()
{
  std::lock_guard<std::mutex> lock(log_lock);
  file = nullptr;
}


std::ostream &
LogStream::get_console()
{
  return *std_out;
}



std::ostringstream &
LogStream::get_stream()
{
  // see if we have already created this stream. if not, do so and
  // set the default flags (why we set these flags is lost to
  // history, but this is what we need to keep several hundred tests
  // from producing different output)
  //
  // note that in all of this we need not worry about thread-safety
  // because we operate on a thread-local object and by definition
  // there can only be one access at a time
  if (outstreams.get().get() == nullptr)
    {
      outstreams.get() = std::make_shared<std::ostringstream>();
      outstreams.get()->setf(std::ios::showpoint | std::ios::left);
    }

  // then return the stream
  return *outstreams.get();
}



std::ostream &
LogStream::get_file_stream()
{
  Assert(file,
         ExcMessage("You can't ask for the std::ostream object for the output "
                    "file if none had been set before."));
  return *file;
}



bool
LogStream::has_file() const
{
  return (file != nullptr);
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
LogStream::push(const std::string &text)
{
  std::string pre;
  if (get_prefixes().size() > 0)
    pre = get_prefixes().top();

  pre += text;
  pre += ":";
  get_prefixes().push(pre);
}



void
LogStream::pop()
{
  if (get_prefixes().size() > 0)
    get_prefixes().pop();
}



std::ios::fmtflags
LogStream::flags(const std::ios::fmtflags f)
{
  return get_stream().flags(f);
}



std::streamsize
LogStream::precision(const std::streamsize prec)
{
  return get_stream().precision(prec);
}



std::streamsize
LogStream::width(const std::streamsize wide)
{
  return get_stream().width(wide);
}



unsigned int
LogStream::depth_console(const unsigned int n)
{
  std::lock_guard<std::mutex> lock(log_lock);
  const unsigned int          h = std_depth;
  std_depth                     = n;
  return h;
}



unsigned int
LogStream::depth_file(const unsigned int n)
{
  std::lock_guard<std::mutex> lock(log_lock);
  const unsigned int          h = file_depth;
  file_depth                    = n;
  return h;
}



bool
LogStream::log_thread_id(const bool flag)
{
  std::lock_guard<std::mutex> lock(log_lock);
  const bool                  h = print_thread_id;
  print_thread_id               = flag;
  return h;
}



std::stack<std::string> &
LogStream::get_prefixes() const
{
  bool                     exists         = false;
  std::stack<std::string> &local_prefixes = prefixes.get(exists);

  // If this is a new locally stored stack, copy the "blessed" prefixes
  // from the initial thread that created logstream.
  if (exists == false)
    {
      const auto parent_prefixes = prefixes.get_for_thread(parent_thread);
      if (parent_prefixes)
        local_prefixes = parent_prefixes.value();
    }

  return local_prefixes;
}



void
LogStream::print_line_head()
{
  const std::string    &head   = get_prefix();
  const std::thread::id thread = std::this_thread::get_id();

  if (get_prefixes().size() <= std_depth)
    {
      if (print_thread_id)
        *std_out << '[' << thread << ']';

      if (head.size() > 0)
        *std_out << head << ':';
    }

  if ((file != nullptr) && (get_prefixes().size() <= file_depth))
    {
      if (print_thread_id)
        *file << '[' << thread << ']';

      if (head.size() > 0)
        *file << head << ':';
    }
}

DEAL_II_NAMESPACE_CLOSE
