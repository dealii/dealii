//----------------------------  logstream.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  logstream.h  ---------------------------
#ifndef __deal2__logstream_h
#define __deal2__logstream_h


#include <base/config.h>
#include <base/exceptions.h>

#include <string>
#include <stack>

// we only need output streams, but older compilers did not provide
// them in a separate include file
#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif


/**
 * A class that simplifies the process of execution logging. It does so by
 * providing
 * @begin{itemize}
 * @item a push and pop mechanism for prefixes, and
 * @item the possibility of distributing information to files and the
 *   console.
 * @end{itemize}
 *
 * The usual usage of this class is through the pregenerated object
 * <tt>deallog</tt>. Typical steps are
 * @begin{itemize}
 * @item <tt>deallog.attach(std::ostream)</tt>: write logging information into a file.
 * @item <tt>deallog.depth_console(n)</tt>: restrict output on screen to outer loops.
 * @item Before entering a new phase of your program, e.g. a new loop,
 *       <tt>deallog.push("loopname")</tt>.
 * @item <tt>deallog << anything << std::endl;</tt> to write logging information
 *       (Usage of <tt>std::endl</tt> is mandatory!).
 * @item <tt>deallog.pop()</tt> when leaving that stage entered with <tt>push</tt>.
 * @end{itemize}
 *
 * @author Guido Kanschat, Wolfgang Bangerth, 1999, 2003
 */
class LogStream
{
  public:
				     /**
				      * Standard constructor, since we
				      * intend to provide an object
				      * <tt>deallog</tt> in the library. Set the
				      * standard output stream to <tt>std::cerr</tt>.
				      */
    LogStream ();

				     /**
				      * Destructor.
				      */
    ~LogStream();
    
				     /**
				      * Enable output to a second
				      * stream <tt>o</tt>.
				      */
    void attach (std::ostream& o);
    
				     /**
				      * Disable output to the second
				      * stream. You may want to call
				      * <tt>close</tt> on the stream that was
				      * previously attached to this object.
				      */
    void detach ();

				     /**
				      * Gives the default stream (<tt>std_out</tt>).
				      */
    std::ostream& get_console ();

				     /**
				      * Gives the file stream.
				      */
    std::ostream& get_file_stream ();

				     /**
				      * Reroutes cerr to LogStream.
				      * Works as a switch, turning
				      * logging of <tt>cerr</tt> on
				      * and off alternatingly with
				      * every call.
				      */
    void log_cerr ();
    
				     /**
				      * Return the prefix string.
				      */
    const std::string& get_prefix () const;
    
				     /**
				      * Push another prefix on the
				      * stack. Prefixes are
				      * automatically separated by a
				      * colon and there is a double
				      * colon after the last prefix.
				      */
    void push (const std::string& text);
        
				     /**
				      * Remove the last prefix.
				      */
    void pop ();
     
				     /**
				      * Maximum number of levels to be
				      * printed on the console. This
				      * function allows to restrict
				      * console output to the upmost
				      * levels of iterations. Only
				      * output with less than <tt>n</tt>
				      * prefixes is printed. By calling
				      * this function with <tt>n=0</tt>, no
				      * console output will be written.
				      *
				      * The previous value of this
				      * parameter is returned.
				      */
    unsigned int depth_console (const unsigned int n);
     
				     /**
				      * Maximum number of levels to be
				      * written to the log file. The
				      * functionality is the same as
				      * <tt>depth_console</tt>, nevertheless,
				      * this function should be used
				      * with care, since it may spoile
				      * the value of a log file.
				      *
				      * The previous value of this
				      * parameter is returned.
				      */
    unsigned int depth_file (const unsigned int n);

                                     /**
				      * Set time printing flag. If this flag
				      * is true, each output line will
				      * be prepended by the user time used
				      * by the running program so far.
				      *
				      * The previous value of this
				      * parameter is returned.
				      */
    bool log_execution_time (const bool flag);

				     /**
				      * Output time differences
				      * between consecutive logs. If
				      * this function is invoked with
				      * <tt>true</tt>, the time difference
				      * between the previous log line
				      * and the recent one is
				      * printed. If it is invoked with
				      * <tt>false</tt>, the accumulated
				      * time since start of the
				      * program is printed (default
				      * behavior).
				      *
				      * The measurement of times is
				      * not changed by this function,
				      * just the output.
				      *
				      * The previous value of this
				      * parameter is returned.
				      */
    bool log_time_differences (const bool flag);
    
				     /**
				      * Output a constant something
				      * through this stream.
				      */
    template <typename T>
    LogStream & operator << (const T &t);

				     /**
				      * Treat ostream
				      * manipulators. This passes on
				      * the whole thing to the
				      * template function with the
				      * exception of the <tt>std::endl</tt>
				      * manipulator, for which special
				      * action is performed: set the
				      * <tt>was_endl</tt> variable that
				      * forces this object to generate
				      * a line head the next time
				      * something is written by this
				      * stream.
				      */
    LogStream & operator<< (std::ostream& (*p) (std::ostream&));

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object. Since sometimes
				      * the size of objects can
				      * not be determined exactly
				      * (for example: what is the
				      * memory consumption of an
				      * STL <tt>std::map</tt> type with a
				      * certain number of
				      * elements?), this is only
				      * an estimate. however often
				      * quite close to the true
				      * value.
				      */
    unsigned int memory_consumption () const;

    				     /**
				      * Exception
				      */
    DeclException0 (ExcNoFileStreamGiven);

  private:
    
				     /**
				      * Stack of strings which are printed
				      * at the beginning of each line to
				      * allow identification where the
				      * output was generated.
				      */
    std::stack<std::string> prefixes;

				     /**
				      * Default stream, where the output
				      * is to go to. This stream defaults
				      * to <tt>std::cerr</tt>, but can be set to another
				      * stream through the constructor.
				      */
    std::ostream  *std_out;

				     /**
				      * Pointer to a stream, where a copy of
				      * the output is to go to. Usually, this
				      * will be a file stream.
				      *
				      * You can set and reset this stream
				      * by the <tt>attach</tt> function.
				      */
    std::ostream  *file;

				     /**
				      * Flag which stores whether the
				      * last operation was a
				      * newline-generation. We use this flag
				      * to generate the list of prefixes at
				      * the next output, rather than
				      * immediately after the newline, since
				      * this might disturb the screen lay-out.
				      */
    bool was_endl;

				     /**
				      * Value denoting the number of
				      * prefixes to be printed to the
				      * standard output. If more than
				      * this number of prefixes is
				      * pushed to the stack, then no
				      * output will be generated until
				      * the number of prefixes shrinks
				      * back below this number.
				      */
    unsigned int std_depth;

				     /**
				      * Same for the maximum depth of
				      * prefixes for output to a file.
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
				      * Original buffer of
				      * <tt>std::cerr</tt>. We store
				      * the address of that buffer
				      * when #log_cerr is called, and
				      * reset it to this value if
				      * #log_cerr is called a second
				      * time, or when the destructor
				      * of this class is run.
				      */
    std::ostream::__streambuf_type* old_cerr;
      
                                     /**
				      * Print head of line. This prints
				      * optional time information and
				      * the contents of the prefix stack.
				      */
    void print_line_head ();

				     /**
				      * Actually do the work of
				      * writing output. This function
				      * unifies the work that is
				      * common to the two
				      * <tt>operator<<</tt> functions.
				      */
    template <typename T>
    void print (const T &t);
};


/* ----------------------------- Inline functions and templates ---------------- */


template <class T>
inline
LogStream &
LogStream::operator<< (const T& t)
{
				   // do the work that is common to
				   // the two operator<< functions
  print (t);
  return *this;
}



inline
LogStream &
LogStream::operator<< (std::ostream& (*p) (std::ostream&))
{
				   // do the work that is common to
				   // the two operator<< functions
  print (p);

				   // next check whether this is the
				   // <tt>endl</tt> manipulator, and if so
				   // set a flag
  std::ostream & (* const p_endl) (std::ostream&) = &std::endl;
  if (p == p_endl)
    was_endl = true;

  return *this;
}



template <class T>
inline
void
LogStream::print (const T &t)
{
				   // if the previous command was an
				   // <tt>std::endl</tt>, print the topmost
				   // prefix and a colon
  if (was_endl)
    {
      print_line_head();
      was_endl = false;
    };

				   // print the rest of the message
  if (prefixes.size() <= std_depth)
    *std_out << t;

  if (file && (prefixes.size() <= file_depth))
    *file << t;
}



/**
 * The standard log object of DEAL.
 *
 * @author Guido Kanschat, 1999
 */
extern LogStream deallog;


#endif
