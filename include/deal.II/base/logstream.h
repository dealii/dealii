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

#ifndef __deal2__logstream_h
#define __deal2__logstream_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/base/thread_local_storage.h>

#include <string>
#include <stack>
#include <map>
#include <cmath>
#include <sstream>

#ifdef HAVE_SYS_TIMES_H
#  include <sys/times.h>
#else
struct tms
{
  int tms_utime, tms_stime, tms_cutime, tms_cstime;
};
#endif


DEAL_II_NAMESPACE_OPEN

/**
 * A class that simplifies the process of execution logging. It does so by
 * providing
 * <ul>
 * <li> a push and pop mechanism for prefixes, and
 * <li> the possibility of distributing information to files and the
 *   console.
 * </ul>
 *
 * The usual usage of this class is through the pregenerated object
 * <tt>deallog</tt>. Typical setup steps are:
 * <ul>
 * <li> <tt>deallog.depth_console(n)</tt>: restrict output on screen to outer loops.
 * <li> <tt>deallog.attach(std::ostream)</tt>: write logging information into a file.
 * <li> <tt>deallog.depth_file(n)</tt>: restrict output to file to outer loops.
 * </ul>
 *
 * Before entering a new phase of your program, e.g. a new loop,
 * a new prefix can be set via <tt>LogStream::Prefix p("loopname");</tt>.
 * The destructor of the prefix will pop the prefix text from the stack.
 *
 * Writes via the <tt>&lt;&lt;</tt> operator,
 * <tt> deallog << "This is a log notice";</tt> will be buffered thread
 * locally until a <tt>std::flush</tt> or <tt>std::endl</tt> is
 * encountered, which will trigger a writeout to the console and, if set
 * up, the log file.
 *
 * <h3>LogStream and thread safety</h3>
 *
 * In the vicinity of concurrent threads, LogStream behaves in the
 * following manner:
 * <ul>
 * <li> Every write to a Logstream with operator <tt>&lt;&lt;</tt> (or with
 * one of the special member functions) is buffered in a thread-local
 * storage.
 * <li> An <tt>std::flush</tt> or <tt>std::endl</tt> will trigger a
 * writeout to the console and (if attached) to the file stream. This
 * writeout is sequentialized so that output from concurrent threads don't
 * interleave.
 * <li> On a new thread, invoking a writeout, as well as a call to #push or
 * #pop will copy the current prefix of the "blessed" thread that created
 * the LogStream instance to a thread-local storage. After that prefixes
 * are thread-local.
 * </ul>
 *
 * <h3>LogStream and reproducible regression test output</h3>
 *
 * Generating reproducible floating point output for regression tests
 * is mildly put a nightmare. In order to make life a little easier,
 * LogStream implements a few features that try to achieve such a
 * goal. These features are turned on by calling test_mode(), and it
 * is not recommended to use them in any other environment. Right now,
 * LogStream implements the following:
 *
 * <ol>
 * <li> A double number very close to zero will end up being output in
 * exponential format, although it has no significant digits. The
 * parameter #double_threshold determines which numbers are too close
 * to zero to be considered nonzero.
 * <li> For float numbers holds the same, but with a typically larger
 * #float_threshold.
 * <li> Rounded numbers become unreliable with inexact
 * arithmetics. Therefore, adding a small number before rounding makes
 * results more reproducible, assuming that numbers like 0.5 are more
 * likely than 0.49997.
 * </ol>
 * It should be pointed out that all of these measures distort the
 * output and make it less accurate. Therefore, they are only
 * recommended if the output needs to be reproducible.
 *
 * @ingroup textoutput
 * @author Guido Kanschat, Wolfgang Bangerth, 1999, 2003, 2011
 */
class LogStream : public Subscriptor
{
public:
  /**
   * A subclass allowing for the safe generation and removal of prefices.
   *
   * Somewhere at the beginning of a block, create one of these objects,
   * and it will appear as a prefix in LogStream output like @p deallog. At
   * the end of the block, the prefix will automatically be removed, when
   * this object is destroyed.
   *
   * In other words, the scope of the object so created determines the
   * lifetime of the prefix. The advantage of using such an object is that
   * the prefix is removed whichever way you exit the scope -- by
   * <code>continue</code>, <code>break</code>, <code>return</code>,
   * <code>throw</code>, or by simply reaching the closing brace. In all of
   * these cases, it is not necessary to remember to pop the prefix
   * manually using LogStream::pop. In this, it works just like the better
   * known Threads::Mutex::ScopedLock class.
   */
  class Prefix
  {
  public:
    /**
     * Set a new prefix for @p deallog, which will be removed when the
     * variable is destroyed.
     */
    Prefix(const std::string &text);

    /**
     * Set a new prefix for the given stream, which will be removed when
     * the variable is destroyed.
     */
    Prefix(const std::string &text,
           LogStream &stream);

    /**
     * Remove the prefix associated with this variable.
     */
    ~Prefix ();

  private:
    SmartPointer<LogStream,LogStream::Prefix> stream;
  };


  /**
   * Standard constructor, since we intend to provide an object
   * <tt>deallog</tt> in the library. Set the standard output stream to
   * <tt>std::cerr</tt>.
   */
  LogStream ();


  /**
   * Destructor.
   */
  ~LogStream();


  /**
   * Enable output to a second stream <tt>o</tt>.
   *
   * The optional argument @p print_job_id specifies whether
   */
  void attach (std::ostream &o,
               const bool    print_job_id = true);


  /**
   * Disable output to the second stream. You may want to call
   * <tt>close</tt> on the stream that was previously attached to this
   * object.
   */
  void detach ();


  /**
   * Setup the logstream for regression test mode.
   *
   * This sets the parameters #double_threshold, #float_threshold, and
   * #offset to nonzero values. The exact values being used have been
   * determined experimentally and can be found in the source code.
   *
   * Called with an argument <tt>false</tt>, switches off test mode and
   * sets all involved parameters to zero.
   */
  void test_mode (bool on=true);


  /**
   * Gives the default stream (<tt>std_out</tt>).
   */
  std::ostream &get_console ();


  /**
   * Gives the file stream.
   */
  std::ostream &get_file_stream ();


  /**
   * @return true, if file stream has already been attached.
   */
  bool has_file () const;


  /**
   * Reroutes cerr to LogStream. Works as a switch, turning logging of
   * <tt>cerr</tt> on and off alternatingly with every call.
   */
  void log_cerr ();


  /**
   * Return the prefix string.
   */
  const std::string &get_prefix () const;


  /**
   * Push another prefix on the stack. Prefixes are automatically separated
   * by a colon and there is a double colon after the last prefix.
   *
   * A simpler way to add a prefix (without the manual need to add the
   * corresponding pop()) is to use the Prefix class.
   */
  void push (const std::string &text);


  /**
   * Remove the last prefix added with push().
   */
  void pop ();


  /**
   * Maximum number of levels to be printed on the console. This function
   * allows to restrict console output to the upmost levels of iterations.
   * Only output with less than <tt>n</tt> prefixes is printed. By calling
   * this function with <tt>n=0</tt>, no console output will be written.
   *
   * The previous value of this parameter is returned.
   */
  unsigned int depth_console (const unsigned int n);


  /**
   * Maximum number of levels to be written to the log file. The
   * functionality is the same as <tt>depth_console</tt>, nevertheless,
   * this function should be used with care, since it may spoile the value
   * of a log file.
   *
   * The previous value of this parameter is returned.
   */
  unsigned int depth_file (const unsigned int n);


  /**
   * Set time printing flag. If this flag is true, each output line will be
   * prepended by the user time used by the running program so far.
   *
   * The previous value of this parameter is returned.
   */
  bool log_execution_time (const bool flag);


  /**
   * Output time differences between consecutive logs. If this function is
   * invoked with <tt>true</tt>, the time difference between the previous
   * log line and the recent one is printed. If it is invoked with
   * <tt>false</tt>, the accumulated time since start of the program is
   * printed (default behavior).
   *
   * The measurement of times is not changed by this function, just the
   * output.
   *
   * The previous value of this parameter is returned.
   */
  bool log_time_differences (const bool flag);


  /**
   * Write detailed timing information.
   */
  void timestamp();


  /**
   * Log the thread id.
   */
  bool log_thread_id (const bool flag);


  /**
   * Set a threshold for the minimal absolute value of double values. All
   * numbers with a smaller absolute value will be printed as zero.
   *
   * The default value for this threshold is zero, i.e. numbers are printed
   * according to their real value.
   *
   * This feature is mostly useful for automated tests: there, one would
   * like to reproduce the exact same solution in each run of a testsuite.
   * However, subtle difference in processor, operating system, or compiler
   * version can lead to differences in the last few digits of numbers, due
   * to different rounding. While one can avoid trouble for most numbers
   * when comparing with stored results by simply limiting the accuracy of
   * output, this does not hold for numbers very close to zero, i.e. zero
   * plus accumulated round-off. For these numbers, already the first digit
   * is tainted by round-off. Using the present function, it is possible to
   * eliminate this source of problems, by simply writing zero to the
   * output in this case.
   */
  void threshold_double(const double t);


  /**
   * The same as threshold_double(), but for float values.
   */
  void threshold_float(const float t);


  /**
   * set the precision for the underlying stream and returns the
   * previous stream precision. This fuction mimics
   * http://www.cplusplus.com/reference/ios/ios_base/precision/
   */
  std::streamsize precision (const std::streamsize prec);


  /**
   * set the width for the underlying stream and returns the
   * previous stream width. This fuction mimics
   * http://www.cplusplus.com/reference/ios/ios_base/width/
   */
  std::streamsize width (const std::streamsize wide);


  /**
   * set the flags for the underlying stream and returns the
   * previous stream flags. This fuction mimics
   * http://www.cplusplus.com/reference/ios/ios_base/flags/
   */
  std::ios::fmtflags flags(const std::ios::fmtflags f);


  /**
   * Output double precision numbers through this stream.
   *
   * If they are set, this function applies the methods for making floating
   * point output reproducible as discussed in the introduction.
   */
  LogStream &operator << (const double t);


  /**
   * Output single precision numbers through this stream.
   *
   * If they are set, this function applies the methods for making floating
   * point output reproducible as discussed in the introduction.
   */
  LogStream &operator << (const float t);


  /**
   * Treat ostream manipulators. This passes on the whole thing to the
   * template function with the exception of the <tt>std::endl</tt>
   * manipulator, for which special action is performed: write the
   * temporary stream buffer including a header to the file and
   * <tt>std::cout</tt> and empty the buffer.
   *
   * An overload of this function is needed anyway, since the compiler
   * can't bind manipulators like @p std::endl directly to template
   * arguments @p T like in the previous general template. This is due to
   * the fact that @p std::endl is actually an overloaded set of functions
   * for @p std::ostream, @p std::wostream, and potentially more of this
   * kind. This function is therefore necessary to pick one element from
   * this overload set.
   */
  LogStream &operator<< (std::ostream& (*p) (std::ostream &));


  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object. Since sometimes the size of objects can not be determined
   * exactly (for example: what is the memory consumption of an STL
   * <tt>std::map</tt> type with a certain number of elements?), this is
   * only an estimate. however often quite close to the true value.
   */
  std::size_t memory_consumption () const;


  /**
   * Exception.
   */
  DeclException0(ExcNoFileStreamGiven);

private:


  /**
   * Internal wrapper around thread-local prefixes. This private
   * function will return the correct internal prefix stack. More
   * important, a new thread-local stack will be copied from the current
   * stack of the "blessed" thread that created this LogStream instance
   * (usually, in the case of deallog, the "main" thread).
   */
  std::stack<std::string> &get_prefixes() const;

  /**
   * Stack of strings which are printed at the beginning of each line to
   * allow identification where the output was generated.
   */
  mutable Threads::ThreadLocalStorage<std::stack<std::string> > prefixes;

  /**
   * Default stream, where the output is to go to. This stream defaults to
   * <tt>std::cerr</tt>, but can be set to another stream through the
   * constructor.
   */
  std::ostream  *std_out;

  /**
   * Pointer to a stream, where a copy of the output is to go to. Usually,
   * this will be a file stream.
   *
   * You can set and reset this stream by the <tt>attach</tt> function.
   */
  std::ostream  *file;

  /**
   * Value denoting the number of prefixes to be printed to the standard
   * output. If more than this number of prefixes is pushed to the stack,
   * then no output will be generated until the number of prefixes shrinks
   * back below this number.
   */
  unsigned int std_depth;

  /**
   * Same for the maximum depth of prefixes for output to a file.
   */
  unsigned int file_depth;

  /**
   * Flag for printing execution time.
   */
  bool print_utime;

  /**
   * Flag for printing time differences.
   */
  bool diff_utime;

  /**
   * Time of last output line.
   */
  double last_time;

  /**
   * Threshold for printing double values. Every number with absolute value
   * less than this is printed as zero.
   */
  double double_threshold;

  /**
   * Threshold for printing float values. Every number with absolute value
   * less than this is printed as zero.
   */
  float float_threshold;

  /**
   * An offset added to every float or double number upon output. This is
   * done after the number is compared to #double_threshold or
   * #float_threshold, but before rounding.
   *
   * This functionality was introduced to produce more reproducible
   * floating point output for regression tests. The rationale is, that an
   * exact output value is much more likely to be 1/8 than 0.124997. If we
   * round to two digits though, 1/8 becomes unreliably either .12 or .13
   * due to machine accuracy. On the other hand, if we add a something
   * above machine accuracy first, we will always get .13.
   *
   * It is safe to leave this value equal to zero. For regression tests,
   * the function test_mode() sets it to a reasonable value.
   *
   * The offset is relative to the magnitude of the number.
   */
  double offset;

  /**
   * Flag for printing thread id.
   */
  bool print_thread_id;

  /**
   * The value times() returned on initialization.
   */
  double reference_time_val;

  /**
   * The tms structure times() filled on initialization.
   */
  struct tms reference_tms;

  /**
   * Original buffer of <tt>std::cerr</tt>. We store the address of that
   * buffer when #log_cerr is called, and reset it to this value if
   * #log_cerr is called a second time, or when the destructor of this
   * class is run.
   */
  std::streambuf *old_cerr;

  /**
   * A flag indicating whether output is currently at a new line
   */
  bool at_newline;

  /**
   * Print head of line. This prints optional time information and the
   * contents of the prefix stack.
   */
  void print_line_head ();

  /**
   * Internal wrapper around "thread local" outstreams. This private
   * function will return the correct internal ostringstream buffer for
   * operater<<.
   */
  std::ostringstream &get_stream();

  /**
   * We use tbb's thread local storage facility to generate a stringstream
   * for every thread that sends log messages.
   */
  Threads::ThreadLocalStorage<std_cxx11::shared_ptr<std::ostringstream> > outstreams;

  template <typename T> friend LogStream &operator << (LogStream &log, const T &t);
};


/* ----------------------------- Inline functions and templates ---------------- */


/**
 * Output a constant something through LogStream:
 *
 * @note We declare this operator as a non-member function so that it is
 * possible to overload it with more specialized templated versions under
 * C++11 overload resolution rules
 */
template <typename T>
inline
LogStream &operator<< (LogStream &log, const T &t)
{
  // print to the internal stringstream
  log.get_stream() << t;
  return log;
}


inline
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
  if (outstreams.get().get() == 0)
    {
      outstreams.get().reset (new std::ostringstream);
      outstreams.get()->setf(std::ios::showpoint | std::ios::left);
    }

  // then return the stream
  return *outstreams.get();
}




inline
LogStream &
LogStream::operator<< (const double t)
{
  std::ostringstream &stream = get_stream();

  // we have to make sure that we don't catch NaN's and +-Inf's with the
  // test, because for these denormals all comparisons are always false.
  // thus, for a NaN, both t<=0 and t>=0 are false at the same time, which
  // can't be said for any other number
  if (! (t<=0) && !(t>=0))
    stream << t;
  else if (std::fabs(t) < double_threshold)
    stream << '0';
  else
    stream << t*(1.+offset);

  return *this;
}



inline
LogStream &
LogStream::operator<< (const float t)
{
  std::ostringstream &stream = get_stream();

  // we have to make sure that we don't catch NaN's and +-Inf's with the
  // test, because for these denormals all comparisons are always false.
  // thus, for a NaN, both t<=0 and t>=0 are false at the same time, which
  // can't be said for any other number
  if (! (t<=0) && !(t>=0))
    stream << t;
  else if (std::fabs(t) < float_threshold)
    stream << '0';
  else
    stream << t*(1.+offset);

  return *this;
}


inline
LogStream::Prefix::Prefix(const std::string &text, LogStream &s)
  :
  stream(&s)
{
  stream->push(text);
}


inline
LogStream::Prefix::~Prefix()
{
  stream->pop();
}


/**
 * The standard log object of deal.II:
 *
 * @author Guido Kanschat, 1999
 */
extern LogStream deallog;


inline
LogStream::Prefix::Prefix(const std::string &text)
  :
  stream(&deallog)
{
  stream->push(text);
}


DEAL_II_NAMESPACE_CLOSE

#endif
