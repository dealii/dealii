// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_logstream_h
#define dealii_logstream_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/thread_local_storage.h>

#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <stack>
#include <string>


DEAL_II_NAMESPACE_OPEN

/**
 * A class that simplifies the process of execution logging. It does so by
 * providing
 * <ul>
 * <li> a push and pop mechanism for prefixes, and
 * <li> the possibility of distributing information to files and the console.
 * </ul>
 *
 * The usual usage of this class is through the pregenerated object
 * <tt>deallog</tt>. Typical setup steps are:
 * <ul>
 * <li> <tt>deallog.depth_console(n)</tt>: restrict output on screen to outer
 * loops.
 * <li> <tt>deallog.attach(std::ostream)</tt>: write logging information into
 * a file.
 * <li> <tt>deallog.depth_file(n)</tt>: restrict output to file to outer
 * loops.
 * </ul>
 *
 * Before entering a new phase of your program, e.g. a new loop, a new prefix
 * can be set via <tt>LogStream::Prefix p("loopname");</tt>. The destructor of
 * the prefix will pop the prefix text from the stack.
 *
 * Write via the <tt>&lt;&lt;</tt> operator, <tt> deallog << "This is a log
 * notice";</tt> will be buffered thread locally until a <tt>std::flush</tt>
 * or <tt>std::endl</tt> is encountered, which will trigger a writeout to the
 * console and, if set up, the log file.
 *
 * <h3>LogStream and thread safety</h3>
 *
 * In the vicinity of concurrent threads, LogStream behaves in the following
 * manner:
 * <ul>
 * <li> Every write to a Logstream with operator <tt>&lt;&lt;</tt> (or with
 * one of the special member functions) is buffered in a thread-local storage.
 * <li> An <tt>std::flush</tt> or <tt>std::endl</tt> will trigger a writeout
 * to the console and (if attached) to the file stream. This writeout is
 * sequentialized so that output from concurrent threads don't interleave.
 * <li> On a new thread, invoking a writeout, as well as a call to #push or
 * #pop will copy the current prefix of the "blessed" thread that created the
 * LogStream instance to a thread-local storage. After that prefixes are
 * thread-local.
 * </ul>
 *
 * @ingroup textoutput
 * @author Guido Kanschat, Wolfgang Bangerth, 1999, 2003, 2011
 */
class LogStream : public Subscriptor
{
public:
  /**
   * A subclass allowing for the safe generation and removal of prefixes.
   *
   * Somewhere at the beginning of a block, create one of these objects, and
   * it will appear as a prefix in LogStream output like @p deallog. At the
   * end of the block, the prefix will automatically be removed, when this
   * object is destroyed.
   *
   * In other words, the scope of the object so created determines the
   * lifetime of the prefix. The advantage of using such an object is that the
   * prefix is removed whichever way you exit the scope -- by
   * <code>continue</code>, <code>break</code>, <code>return</code>,
   * <code>throw</code>, or by simply reaching the closing brace. In all of
   * these cases, it is not necessary to remember to pop the prefix manually
   * using LogStream::pop(). In this, it works just like the better known
   * std::unique_ptr and std::lock_guard classes.
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
     * Set a new prefix for the given stream, which will be removed when the
     * variable is destroyed.
     */
    Prefix(const std::string &text, LogStream &stream);

    /**
     * Remove the prefix associated with this variable.
     */
    ~Prefix();

  private:
    /**
     * A pointer to the LogStream object to which the prefix is
     * applied.
     */
    SmartPointer<LogStream, LogStream::Prefix> stream;
  };


  /**
   * Standard constructor. The standard output stream to
   * <tt>std::cout</tt>.
   */
  LogStream();


  /**
   * Destructor.
   */
  ~LogStream() override;


  /**
   * Enable output to a second stream <tt>o</tt>.
   *
   * @param o Attach this output stream.
   *
   * @param[in] print_job_id Whether or not the JobIdentifier for the current
   * process should be printed to the stream.
   *
   * @param[in] flags Format flags to set on the output stream @p o.
   */
  void
  attach(std::ostream &                o,
         const bool                    print_job_id = true,
         const std::ios_base::fmtflags flags        = std::ios::showpoint |
                                               std::ios::left);


  /**
   * Disable output to the second stream. You may want to call <tt>close</tt>
   * on the stream that was previously attached to this object.
   */
  void
  detach();


  /**
   * Return the default stream (<tt>std_out</tt>).
   */
  std::ostream &
  get_console();


  /**
   * Return the file stream.
   */
  std::ostream &
  get_file_stream();


  /**
   * Return @p true if file stream has already been attached,
   * @p false otherwise.
   */
  bool
  has_file() const;


  /**
   * Return the prefix string.
   */
  const std::string &
  get_prefix() const;


  /**
   * Push another prefix on the stack. Prefixes are automatically separated by
   * a colon and there is a double colon after the last prefix.
   *
   * A simpler way to add a prefix (without the manual need to add the
   * corresponding pop()) is to use the LogStream::Prefix class. Using
   * that class has the advantage that the corresponding pop() call is
   * issued whenever the Prefix object goes out of scope -- either at
   * the end of the code block, at the nearest @p return statement, or
   * because an intermediate function call results in an exception that
   * is not immediately caught.
   */
  void
  push(const std::string &text);


  /**
   * Remove the last prefix added with push().
   */
  void
  pop();


  /**
   * Maximum number of levels to be printed on the console. The default is 0,
   * which will not generate any output. This function allows one to restrict
   * console output to the highest levels of iterations. Only output with less
   * than <tt>n</tt> prefixes is printed. By calling this function with
   * <tt>n=0</tt>, no console output will be written. See step-3 for an
   * example usage of this method.
   *
   * The previous value of this parameter is returned.
   */
  unsigned int
  depth_console(const unsigned int n);


  /**
   * Maximum number of levels to be written to the log file. The functionality
   * is the same as <tt>depth_console</tt>, nevertheless, this function should
   * be used with care, since it may spoil the value of a log file.
   *
   * The previous value of this parameter is returned.
   */
  unsigned int
  depth_file(const unsigned int n);


  /**
   * Log the thread id.
   */
  bool
  log_thread_id(const bool flag);


  /**
   * set the precision for the underlying stream and returns the previous
   * stream precision. This function mimics
   * http://www.cplusplus.com/reference/ios/ios_base/precision/
   */
  std::streamsize
  precision(const std::streamsize prec);


  /**
   * set the width for the underlying stream and returns the previous stream
   * width. This function mimics
   * http://www.cplusplus.com/reference/ios/ios_base/width/
   */
  std::streamsize
  width(const std::streamsize wide);


  /**
   * set the flags for the underlying stream and returns the previous stream
   * flags. This function mimics
   * http://www.cplusplus.com/reference/ios/ios_base/flags/
   */
  std::ios::fmtflags
  flags(const std::ios::fmtflags f);


  /**
   * Treat ostream manipulators. This passes on the whole thing to the
   * template function with the exception of the <tt>std::endl</tt>
   * manipulator, for which special action is performed: write the temporary
   * stream buffer including a header to the file and <tt>std::cout</tt> and
   * empty the buffer.
   *
   * An overload of this function is needed anyway, since the compiler can't
   * bind manipulators like @p std::endl directly to template arguments @p T
   * like in the previous general template. This is due to the fact that @p
   * std::endl is actually an overloaded set of functions for @p std::ostream,
   * @p std::wostream, and potentially more of this kind. This function is
   * therefore necessary to pick one element from this overload set.
   */
  LogStream &
  operator<<(std::ostream &(*p)(std::ostream &));


  /**
   * Return an estimate for the memory consumption, in bytes, of this object.
   * This is not exact (but will usually be close) because calculating the
   * memory usage of trees (e.g., <tt>std::map</tt>) is difficult.
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * Internal wrapper around thread-local prefixes. This private function will
   * return the correct internal prefix stack. More important, a new thread-
   * local stack will be copied from the current stack of the "blessed" thread
   * that created this LogStream instance (usually, in the case of deallog,
   * the "main" thread).
   */
  std::stack<std::string> &
  get_prefixes() const;

  /**
   * Stack of strings which are printed at the beginning of each line to allow
   * identification where the output was generated.
   */
  mutable Threads::ThreadLocalStorage<std::stack<std::string>> prefixes;

  /**
   * We record the thread id of the thread creating this object. We need
   * this information to "steal" the current prefix from this "master"
   * thread on first use of deallog on a new thread.
   */
  std::thread::id master_thread;

  /**
   * Default stream, where the output is to go to. This stream defaults to
   * <tt>std::cout</tt>, but can be set to another stream through the
   * constructor.
   */
  std::ostream *std_out;

  /**
   * Pointer to a stream, where a copy of the output is to go to. Usually,
   * this will be a file stream.
   *
   * You can set and reset this stream by the <tt>attach</tt> function.
   */
  std::ostream *file;

  /**
   * Value denoting the number of prefixes to be printed to the standard
   * output. If more than this number of prefixes is pushed to the stack, then
   * no output will be generated until the number of prefixes shrinks back
   * below this number.
   */
  unsigned int std_depth;

  /**
   * Same for the maximum depth of prefixes for output to a file.
   */
  unsigned int file_depth;

  /**
   * Flag for printing thread id.
   */
  bool print_thread_id;

  /**
   * A flag indicating whether output is currently at a new line
   */
  bool at_newline;

  /**
   * Print head of line.
   */
  void
  print_line_head();

  /**
   * Internal wrapper around "thread local" outstreams. This private function
   * will return the correct internal ostringstream buffer for operator<<.
   */
  std::ostringstream &
  get_stream();

  /**
   * We use our thread local storage facility to generate a stringstream for
   * every thread that sends log messages.
   */
  Threads::ThreadLocalStorage<std::shared_ptr<std::ostringstream>> outstreams;

  template <typename T>
  friend LogStream &
  operator<<(LogStream &log, const T &t);
};


/* ----------------------------- Inline functions and templates ----------------
 */


/**
 * Output a constant something through LogStream:
 *
 * @note We declare this operator as a non-member function so that it is
 * possible to overload it with more specialized templated versions under
 * C++11 overload resolution rules
 */
template <typename T>
inline LogStream &
operator<<(LogStream &log, const T &t)
{
  // print to the internal stringstream
  log.get_stream() << t;
  return log;
}



/**
 * The standard log object of deal.II:
 *
 * @author Guido Kanschat, 1999
 */
extern LogStream deallog;



DEAL_II_NAMESPACE_CLOSE

#endif
