//----------------------------  logstream.h  ---------------------------
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
//----------------------------  logstream.h  ---------------------------
#ifndef __deal2__logstream_h
#define __deal2__logstream_h


#include <base/exceptions.h>

#include <string>
#include <stack>
#include <iostream>


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
 * @p{deallog}. Typical steps are
 * <OL>
 * <LI> @p{deallog.attach(ostream)}: write logging information into a file.
 * <LI> @p{deallog.depth_console(n)}: restrict output on screen to outer loops.
 * <LI> Before entering a new phase of your program, e.g. a new loop,
 * @p{deallog.push("loopname")}.
 * <LI> @p{deallog << anything << endl;} to write logging information
 * (Usage of @p{endl} is mandatory!).
 * <LI> @p{deallog.pop()} when leaving that stage entered with @p{push}.
 * </OL>
 *
 * @author Guido Kanschat, Wolfgang Bangerth, 1999
 */
class LogStream
{
  private:
    
				     /**
				      * Stack of strings which are printed
				      * at the beginning of each line to
				      * allow identification where the
				      * output was generated.
				      */
    stack<string> prefixes;

				     /**
				      * Default stream, where the output
				      * is to go to. This stream defaults
				      * to @p{cerr}, but can be set to another
				      * stream through the constructor.
				      */
    ostream  *std_out;

				     /**
				      * Pointer to a stream, where a copy of
				      * the output is to go to. Usually, this
				      * will be a file stream.
				      *
				      * You can set and reset this stream
				      * by the @p{attach} function.
				      */
    ostream  *file;

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

  public:
				     /**
				      * Standard constructor, since we
				      * intend to provide an object
				      * @p{deallog} in the library. Set the
				      * standard output stream to @p{cerr}.
				      */
    LogStream();
    
				     /**
				      * Enable output to a second
				      * stream @p{o}.
				      */
    void attach(ostream& o);
    
				     /**
				      * Disable output to the second
				      * stream. You may want to call
				      * @p{close} on the stream that was
				      * previously attached to this object.
				      */
    void detach();

				     /**
				      * Gives the default stream (@p{std_out}).
				      */
    ostream& get_console();

				     /**
				      * Gives the file stream.
				      */
    ostream& get_file_stream();
    
				     /**
				      * Push another prefix on the
				      * stack. Prefixes are
				      * automatically separated by a
				      * colon and there is a double
				      * colon after the last prefix.
				      */
    void push (const string& text);
        
				     /**
				      * Remove the last prefix.
				      */
    void pop();
     
				     /**
				      * Maximum number of levels to be
				      * printed on the console. This
				      * function allows to restrict
				      * console output to the upmost
				      * levels of iterations. Only
				      * output with less than @p{n}
				      * prefixes is printed. By calling
				      * this function with @p{n=0}, no
				      * console output will be written.
				      */
    void depth_console (unsigned int n);
     
				     /**
				      * Maximum number of levels to be
				      * written to the log file. The
				      * functionality is the same as
				      * @p{depth_console}, nevertheless,
				      * this function should be used
				      * with care, since it may spoile
				      * the value of a log file.
				      */
    void depth_file (unsigned int n);

                                     /**
				      * Set time printing flag. If this flag
				      * is true, each output line will
				      * be prepended by the user time used
				      * by the running program so far.
				      */
    void log_execution_time(bool flag);

				     /**
				      * Output a constant something through
				      * this stream.
				      */
    template <typename T>
    LogStream & operator << (const T &t);

				     /**
				      * Output a function. This really is not
				      * to output the function, but calls the
				      * function on the present object. It
				      * works the same way as it is possible
				      * to output the stream manipulators
				      * (like @p{endl}, etc) work on standard
				      * streams: @p{stream << endl;}. We need
				      * to overload this function in order
				      * to know when we have to print the
				      * prefixes, since a user may use several
				      * calls to @p{operator <<} before the
				      * line is complete. The overloaded
				      * function @p{void endl(LogStream &)}
				      * tells the @p{LogStream} that the end of
				      * the line is reached.
				      */
    LogStream & operator << (void (f)(LogStream &));

				     /**
				      * Declare this function as a friend.
				      */
    friend void endl (LogStream &);

    				     /**
				      * Exception
				      */
    DeclException0 (ExcNoFileStreamGiven);

  private:
                                     /**
				      * Print head of line. This prints
				      * optional time information and
				      * the contents of the prefix stack.
				      */
    void print_line_head();
};


/* ----------------------------- Inline functions and templates ---------------- */


template <class T>
inline
LogStream &
LogStream::operator<< (const T& t)
{
				   // if the previous command was an
				   // @p{endl}, print the topmost prefix
				   // and a colon
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

  return *this;
}


/**
 * Replacement of @p{endl} for @p{LogStream}.
 *
 * Overloaded version of the stream manipulator function @p{endl} which
 * results in calling the original version of @p{endl} for each of the
 * two streams, if the present prefix number does not exceed the
 * specified maximal number.
 *
 * @author Guido Kanschat, 1999
 */
inline void endl(LogStream& s)
{
  if (s.prefixes.size() <= s.std_depth)
    *s.std_out << endl;

  if (s.file && (s.prefixes.size() <= s.file_depth))
    *s.file << endl;

  s.was_endl = true;
};


/**
 * The standard log object of DEAL.
 *
 * @author Guido Kanschat, 1999
 */
extern LogStream deallog;


#endif

