/*----------------------------   logstream.h     ---------------------------*/
/*      $Id$                 */
#ifndef __logstream_H
#define __logstream_H
/*----------------------------   logstream.h     ---------------------------*/


#include <string>
#include <stack>
#include <iostream>




/**
 * A class that simplifies the process of execution logging. It does so by
 * providing
 * \begin{itemize}
 * \item a push and pop mechanism for prefixes, and
 * \item the possibility of distributing information to files and the
 *   console.
 * \end{itemize}
 *
 * The usual usage of this class is through the pregenerated object
 * #deallog#. Typical steps are
 * <OL>
 * <LI> #deallog.attach(ostream)#: write logging information into a file.
 * <LI> #deallog.depth_console(n)#: restrict output on screen to outer loops.
 * <LI> Before entering a new phase of your program, e.g. a new loop,
 * #deallog.push("loopname")#.
 * <LI> #deallog << anything << endl;# to write logging information
 * (Usage of #endl# is mandatory!).
 * <LI> #deallog.pop()# when leaving that stage entered with #push#.
 * </OL>
 *
 * @author Guido Kanschat, 1999
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
				      * to #cerr#, but can be set to another
				      * stream through the constructor.
				      */
    ostream  *std_out;

				     /**
				      * Pointer to a stream, where a copy of
				      * the output is to go to. Usually, this
				      * will be a file stream.
				      *
				      * You can set and reset this stream
				      * by the #attach# function.
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
    
  public:
				     /**
				      * Standard constructor, since we
				      * intend to provide an object
				      * #deallog# in the library. Set the
				      * standard output stream to #cerr#.
				      */
    LogStream();
    
				     /**
				      * Enable output to a second
				      * stream #o#.
				      */
    void attach(ostream& o);
    
				     /**
				      * Disable output to the second
				      * stream. You may want to call
				      * #close# on the stream that was
				      * previously attached to this object.
				      */
    void detach();
    

    
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
				      * output with less than #n#
				      * prefixes is printed. By calling
				      * this function with #n=0#, no
				      * console output will be written.
				      */
    void depth_console (unsigned int n);
     
				     /**
				      * Maximum number of levels to be
				      * written to the log file. The
				      * functionality is the same as
				      * #depth_console#, nevertheless,
				      * this function should be used
				      * with care, since it may spoile
				      * the value of a log file.
				      */
    void depth_file (unsigned int n);

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
				      * (like #endl#, etc) work on standard
				      * streams: #stream << endl;#. We need
				      * to overload this function in order
				      * to know when we have to print the
				      * prefixes, since a user may use several
				      * calls to #operator <<# before the
				      * line is complete. The overloaded
				      * function #void endl(LogStream &)#
				      * tells the #LogStream# that the end of
				      * the line is reached.
				      */
    LogStream & operator << (void (f)(LogStream &));

				     /**
				      * Declare this function as a friend.
				      */
    friend void endl (LogStream &);
};




/* ----------------------------- Inline functions and templates ---------------- */


inline
void
LogStream::push (const string& text)
{
				   // strange enough: if make this
				   // function non-inline with
				   // gcc2.8, we get very strange
				   // compiler errors...
  string pre=prefixes.top();
  pre += text;
  pre += string(":");
  prefixes.push(pre);
}



// sorry for the weird following declaration spanning 5 lines, but
// doc++ gets confused if we be more compact :-(

template <class T>
inline
LogStream &
LogStream::
operator<< (const T& t)
{
				   // if the previous command was an
				   // #endl#, print the topmost prefix
				   // and a colon
  if (was_endl)
    {
      if (prefixes.size() <= std_depth)
	*std_out << prefixes.top() << ':';

      if (file && (prefixes.size() <= file_depth))
	*file << prefixes.top() << ':';

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
 * Replacement of #endl# for #LogStream#.
 *
 * Overloaded version of the stream manipulator function #endl# which
 * results in calling the original version of #endl# for each of the
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




/*----------------------------   logstream.h     ---------------------------*/
/* end of #ifndef __logstream_H */
#endif
/*----------------------------   logstream.h     ---------------------------*/

