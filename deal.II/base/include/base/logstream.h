/*----------------------------   logstream.h     ---------------------------*/
/*      $Id$                 */
#ifndef __logstream_H
#define __logstream_H
/*----------------------------   logstream.h     ---------------------------*/


#include <stack>
#include <string>
#include <iostream>



/**
 * A class that simplifies the process of execution logging by
 * providing
 * \begin{itemize}
 * \item a push and pop mechanism for prefixes and
 * \item the possibility of distributing information to files and the
 *   console.
 * \end{itemize}
 */
class LogStream
{
public:
  stack<string> prefixes;
  ostream& std;
  ostream* file;

  bool was_endl;
  unsigned std_depth, file_depth;
  
public:
				   /**
				    * Standard constructor, since we
				    * intend to provide an object
				    * clog in the library.
				    */
  LogStream();
    
  				   /**
				    * Enable output to a second
				    * stream.
				    */
  void attach(ostream& o);
    
				   /**
				    * Disable output to the second
				    * stream. 
				    */
  void detach()
      {
	file = 0;
      }
    
				     /**
				      * Push another prefix on the
				      * stack. Prefixes are
				      * automatically separated by a
				      * colon and there is a double
				      * colon after the last prefix.
				      */
    void push(const string& text)
      {
	string pre=prefixes.top();
	pre += text;
	pre += string(":");
	prefixes.push(pre);
      }
    
				     /**
				      * Remove the last prefix.
				      */
    void pop()
      {
  	prefixes.pop();
      }
     
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
    void depth_console(unsigned n)
      {
  	std_depth = n;
      }
     
				     /**
				      * Maximum number of levels to be
				      * written to the log file. The
				      * functionality is the same as
				      * #depth_console#, nevertheless,
				      * this function should be used
				      * with care, since it may spoile
				      * the value of a log file.
				      */
    void depth_file(unsigned n)
      {
 	file_depth = n;
      }

    ostream & get_std_stream () const {
      return std;
    };
};




/* ----------------------------- Inline functions and templates ---------------- */


template <class T>
inline void
writestuff(LogStream& s, const T& t)
{
  if (s.prefixes.size()<=s.std_depth)
    s.get_std_stream () << t;
  if (s.file && (s.prefixes.size()<=s.file_depth))
    *(s.file) << t;
}




template <class T>
inline LogStream&
operator << (LogStream& s, const T& t)
{
  if (s.was_endl)
  {
    writestuff(s, s.prefixes.top());
    writestuff(s,':');
  }
  s.was_endl = false;
  writestuff(s,t);
  return s;
}

//template <>
inline LogStream&
operator << (LogStream& s, const char* c)
{
  if (s.was_endl)
  {
    writestuff(s, s.prefixes.top());
    writestuff(s,':');
  }
  s.was_endl = false;
  writestuff(s, c);
  return s;
}

inline void endl(LogStream& s)
{
  ostream& (*p) (ostream &) = &endl;
  writestuff(s, p);
  s.was_endl = true;
}

//template <>
inline LogStream&
operator << (LogStream& s, void (*f)(LogStream&))
{
  f(s);
  return s;
}

extern LogStream deallog;




/*----------------------------   logstream.h     ---------------------------*/
/* end of #ifndef __logstream_H */
#endif
/*----------------------------   logstream.h     ---------------------------*/
