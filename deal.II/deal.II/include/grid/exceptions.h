/*----------------------------   exceptions.h     ---------------------------*/
/*      $Id$  */
#ifndef __exceptions_H
#define __exceptions_H
/*----------------------------   exceptions.h     ---------------------------*/


#include <iostream.h>



#ifndef __GNUC__
#  define __PRETTY_FUNCTION__ "(unknown)"
#endif

				 /**
				  *  This class should be used as a base class for
				  *  all exception classes. Do not use its methods
				  *  and variables directly since the interface
				  *  and mechanism may be subject to change. Rather
				  *  create new exception classes using the
				  *  #DeclException# macro family.
				  *
				  *  @see Assert
				  *  @see DeclException0
				  *  @author Wolfgang Bangerth, November 1997
				  */
class ExceptionBase {
  public:
    ExceptionBase () {};
				     /**
				      *  The constructor takes the file in which the
				      *  error happened, the line and the violated
				      *  condition as well as the name of the
				      *  exception class as a #char*# as arguments.
				      */
    ExceptionBase (const char* f, const int l, const char *func,
		   const char* c, const char *e) :
		    file(f), line(l), function(func), cond(c), exc(e) {};

				     /**
				      *  Set the file name and line of where the
				      *  exception appeared as well as the violated
				      *  condition and the name of the exception as
				      *  a char pointer.
				      */
    void SetFields (const char* f,
		    const int l,
		    const char *func,
		    const char *c,
		    const char *e) {
      file = f;
      line = l;
      function = func;
      cond = c;
      exc  = e;
    };
				     /**
				      *  Print out the general part of the error
				      *  information.
				      */
    void PrintExcData (ostream &out) const {
      out << "An error occurred in line <" << line
	  << "> of file <" << file
	  << "> in function" << endl
	  << "    " << function << endl
	  << "The violated condition was: "<< endl
	  << "    " << cond << endl
	  << "The name and call sequence of the exception was:" << endl
	  << "    " << exc  << endl
	  << "Additional Information: " << endl;
    };

				     /**
				      *  Print more specific information about the
				      *  exception which occured. Overload this
				      *  function in your own exception classes.
				      *  Since we use templates rather than derivation
				      *  information, we need not declare this
				      *  function as #virtual# and thus avoid problems
				      *  with virtual functions in classes in which
				      *  all functions are declared inline.
				      */
    void PrintInfo (ostream &out) const {
      out << "(none)" << endl;
    };

    
    const char *file;
    int         line;
    const char *function;
    const char *cond;
    const char *exc;
};


				 /**
				  *  This routine does the main work for the
				  *  exception generation mechanism.
				  *
				  *  @see Assert
				  */
template <class exc>
void __IssueError (const char *file,
		   int         line,
		   const char *function,
		   const char *cond,
		   const char *exc,
		   exc         e) {
				   // Fill the fields of the exception object
      e.SetFields (file, line, function, cond, exc);
      cerr << "--------------------------------------------------------"
	   << endl;
				       // print out general data
      e.PrintExcData (cerr);
				       // print out exception specific data
      e.PrintInfo (cerr);
      cerr << "--------------------------------------------------------"
	   << endl;
      
      abort ();
    };


    
#ifdef DEBUG  ////////////////////////////////////////

/**
    This is the main routine in the exception mechanism.
    Use it in the following way: declare an exception
    class using the #DeclException# macros (see there),
    for example
    \begin{verbatim}
      DeclException2 (ExcDomain, int, int,
                      << "Index= " << arg1 << "Upper Bound= " << arg2);
    \end{verbatim}
    to declare an exception class named #ExcDomain#, which
    has two variables as additional information (named
    #arg1# and #arg2# by default) and which outputs the
    given sequence (which is appended to an #ostream#
    variable's name, thus the weird syntax). There are
    other #DeclExceptionN# macros for exception classes
    with more or no parameters. It is proposed to let start
    the name of all exceptions by #Exc...#.

    In your source code write lines like
    \begin{verbatim}
      Assert (n<dim, ExcDomain(n,dim));
    \end{verbatim}
    which does the following:
    \begin{verbatim}
      #ifdef DEBUG
          if (!cond)
                issue error of class ExcDomain(n,dim)
      #else
          do nothing
      #endif
    \end{verbatim}
    i.e. it issues an error only if the preprocessor variable
    #DEBUG# is set and if the given condition (in this case
    #n<dim# is violated.

    If the exception was declared using the #DeclException0 (...)#
    macro, i.e. without any additional parameters, its name has
    nonetheless to be given with parentheses:
    #Assert (i>m, ExcSomewhat());#

    {\bf How it works internally:}
    The call #Assert (cond, exc);#
    is converted by the preprocessor into the sequence
    \begin{verbatim}
      if (!(cond))
        __IssueError (__FILE__,
                      __LINE__,
                      __PRETTY_FUNCTION__,
                      #cond,
                      #exc,
                      &exc);
    \end{verbatim}
    i.e. if the given condition is violated, then the file and
    line in which the exception occured as well as
    the condition itself and the call sequence of the
    exception object is transferred. Additionally an object
    of the form given by #exc# is created (this is normally an
    unnamed object like in #ExcDomain (n, dim)# of class
    #ExcDomain#) and transferred to the #__IssueError#
    function.

    #__PRETTY__FUNCTION__# is a macro defined only by the GNU CC
    compiler and gives the name of the function. If another compiler
    is used, we set #__PRETTY_FUNCTION__ = "(unknown)"#.

    In #__IssueError# the given data
    is transferred into the #exc# object by calling the
    #SetFields# function; after that, the general error info
    is printed onto #cerr# using the #PrintError# function of
    #exc# and finally the exception specific data is printed
    using the user defined function #PrintError# (which is
    normally created using the #DeclException (...)# macro
    family.

    After printing all this information, #abort()# is called.
    This may in future versions be replaced by calling
    #throw exc;#.

    If the preprocessor variable #DEBUG# is not set, then nothing
    happens, i.e. the macro is expanded to #(void) 0;#.

    @memo Assert that a certain condition is fulfilled, otherwise
    issue an error and abort the program.
    
    @see DeclException0
    @author Wolfgang Bangerth, November 1997
    */
#define Assert(cond, exc)                      \
  if (!(cond))                                 \
    __IssueError (__FILE__,                    \
		  __LINE__,                    \
		  __PRETTY_FUNCTION__, #cond, #exc, exc)


#else        ////////////////////////////////////////

#define Assert(cond, exc)                     \
  (void) 0
#endif      ////////////////////////////////////////



/**
    Declare an exception class without any additional parameters.
    There is a whole family of #DeclException?# macros
    where #?# is to be replaced by the number of additional
    parameters (0 to 5 presently).
    
    The syntax is as follows:
    \begin{verbatim}
      DeclException2 (ExcDomain,
                      int,
                      int,
                      << " i=" << arg1 << ", m=" << arg2);
    \end{verbatim}
    The first is the name of the exception class to be created.
    The next arguments are the types of the parameters (in this
    case there are two type names needed) and finally the output
    sequence with which you can print additional information.
    
    The syntax of the output sequence is a bit weird but gets
    clear once you see how this macro is defined:
    \begin{verbatim}
    class name : public ExceptionBase {
      public:
        name (const type1 a1, const type2 a2) :
                       arg1 (a1), arg2(a2) {};
        void PrintInfo (ostream &out) const {
          out outsequence << endl;
        };
      private:
        type1 arg1;
        type2 arg2;
    };
    \end{verbatim}
    The #PrintInfo# function is not declared #virtual# since we would
    then create a class with virtual functions and only inlined code.
    For this kind of class at least GNU C++ does not know where to put
    the virtual method table and the program will not be compiled.
    
    If declared as specified, you can later use this exception class
    in the following manner:
    \begin{verbatim}
      int i=5;
      int m=3;
      Assert (i<m, MyExc2(i,m));
    \end{verbatim}
    and the output if the condition fails will be
    \begin{verbatim}
      --------------------------------------------------------
      An error occurred in line <301> of file <exc-test.cc>.
      The violated condition was: 
        i<m
      The name and call sequence of the exception was:
        MyExc2(i,m)
      Additional Information: 
        i=5, m=3
      --------------------------------------------------------
    \end{verbatim}
    
    Obviously for the #DeclException0(name)# macro, no types and
    also no output sequence is allowed.
    
    @author Wolfgang Bangerth, November 1997
    */
#define DeclException0(Exception0)  \
class Exception0 :  public ExceptionBase {};

				 /**
				  *  Declare an exception class with
				  *  one additional parameter.
				  *
				  *  @see DeclException0
				  */
#define DeclException1(Exception1, type1, outsequence)                \
class Exception1 : public ExceptionBase {                             \
  public:                                                             \
      Exception1 (const type1 a1) : arg1 (a1) {};                     \
      void PrintInfo (ostream &out) const {                           \
        out outsequence << endl;                                      \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
}

				 /**
				  *  Declare an exception class with
				  *  two additional parameters.
				  *
				  *  @see DeclException0
				  */
#define DeclException2(Exception2, type1, type2, outsequence)         \
class Exception2 : public ExceptionBase {                             \
  public:                                                             \
      Exception2 (const type1 a1, const type2 a2) :                   \
	      arg1 (a1), arg2(a2) {};                                 \
      void PrintInfo (ostream &out) const {                           \
        out outsequence << endl;                                      \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
      const type2 arg2;                                               \
}

				 /**
				  *  Declare an exception class with
				  *  three additional parameters.
				  *
				  *  @see DeclException0
				  */
#define DeclException3(Exception3, type1, type2, type3, outsequence)  \
class Exception3 : public ExceptionBase {                             \
  public:                                                             \
      Exception3 (const type1 a1, const type2 a2, const type3 a3) :   \
	      arg1 (a1), arg2(a2), arg3(a3) {};                       \
      void PrintInfo (ostream &out) const {                           \
        out outsequence << endl;                                      \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
      const type2 arg2;                                               \
      const type3 arg3;                                               \
}

				 /**
				  *  Declare an exception class with
				  *  four additional parameters.
				  *
				  *  @see DeclException0
				  */
#define DeclException4(Exception4, type1, type2, type3, type4, outsequence) \
class Exception4 : public ExceptionBase {                             \
  public:                                                             \
      Exception4 (const type1 a1, const type2 a2,                     \
	    const type3 a3, const type4 a4) :                         \
	      arg1 (a1), arg2(a2), arg3(a3), arg4(a4) {};             \
      void PrintInfo (ostream &out) const {                           \
        out outsequence << endl;                                      \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
      const type2 arg2;                                               \
      const type3 arg3;                                               \
      const type4 arg4;                                               \
}

				 /**
				  *  Declare an exception class with
				  *  five additional parameters.
				  *
				  *  @see DeclException0
				  */
#define DeclException5(Exception5, type1, type2, type3, type4, type5, outsequence) \
class Exception5 : public ExceptionBase {                             \
  public:                                                             \
      Exception5 (const type1 a1, const type2 a2, const type3 a3,     \
	    const type4 a4, const type5 a5) :                         \
	      arg1 (a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {};   \
      void PrintInfo (ostream &out) const {                           \
        out outsequence << endl;                                      \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
      const type2 arg2;                                               \
      const type3 arg3;                                               \
      const type4 arg4;                                               \
      const type5 arg5;                                               \
}






/*----------------------------   exceptions.h     ---------------------------*/
/* end of #ifndef __exceptions_H */
#endif
/*----------------------------   exceptions.h     ---------------------------*/
