//----------------------------  exceptions.h  ---------------------------
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
//----------------------------  exceptions.h  ---------------------------
#ifndef __deal2__exceptions_h
#define __deal2__exceptions_h


#include <iostream>
#include <exception>
#include <cstdlib>

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
 *
 *  \section{General overview of the exception handling mechanism in #deal.II#}
 *
 *  The error handling mechanism in #deal.II# is generally used in two ways.
 *  The first uses error checking in debug mode only and is useful for programs
 *  which are not fully tested. When the program shows no error anymore, one may
 *  switch off error handling and get better performance by this, since checks
 *  for errors are done quite frequently in the library (a typical speed up is
 *  a factor of four!). This mode of exception generation is most useful for
 *  internal consistency checks such as range checking or checking of the
 *  validity of function arguments. Errors of this kind usually are programming
 *  errors and the program should abort with as detailed a message as possible,
 *  including location and reason for the generation of the exception.
 *
 *  The second mode is for error checks which should always be on, such as for
 *  I/O errors, failing memory requests and the like. It does not make much
 *  sense to turn this mode off, since this kind of errors may happen in tested
 *  and untested programs likewise. Exceptions of this kind do not terminate the
 *  program, rather they throw exceptions in the #C++# manner, allowing the
 *  program to catch them and eventually do something about it. As it may be
 *  useful to have some information printed out if an exception could not be
 *  handled properly, additional information is passed along as for the first
 *  mode. The latter makes it necessary to provide a family of macros which
 *  enter this additional information into the exception class; this could
 *  in principle be done by the programmer himself each time by hand, but since
 *  the information can be obtained automatically, a macro is provided for
 *  this.
 *
 *  Both modes use exception classes, which need to have special features
 *  additionally to the #C++# standard's #exception# class. Such a class
 *  is declared by the following lines of code:
 *   \begin{verbatim}
 *     DeclException2 (ExcDomain, int, int,
 *                     << "Index= " << arg1 << "Upper Bound= " << arg2);
 *  \end{verbatim}
 *  This declares an exception class named #ExcDomain#, which
 *  has two variables as additional information (named
 *  #arg1# and #arg2# by default) and which outputs the
 *  given sequence (which is appended to an #ostream#
 *  variable's name, thus the weird syntax). There are
 *  other #DeclExceptionN# macros for exception classes
 *  with more or no parameters. It is proposed to let start
 *  the name of all exceptions by #Exc...# and to declare them locally to
 *  the class it is to be used in. Declaring exceptions globally is possible
 *  but pollutes the global namespace, is less readable and thus unnecessary.
 *
 *  Since exception classes are declared the same way for both modes of error
 *  checking, it is possible to use an exception declared through the
 *  #DeclExceptionN(...)# macro family in both modes; there is no need to
 *  declare different classes for each of these.
 *
 *
 *  \section{Use of the debug mode exceptions}
 *
 *  To use the exception mechanism for debug mode error checking, write lines
 *  like the following in your source code:
 *  \begin{verbatim}
 *    Assert (n<dim, ExcDomain(n,dim));
 *  \end{verbatim}
 *  which by macro expansion does the following:
 *  \begin{verbatim}
 *    #ifdef DEBUG
 *        if (!cond)
 *              issue error of class ExcDomain(n,dim)
 *    #else
 *        do nothing
 *    #endif
 *  \end{verbatim}
 *  i.e. it issues an error only if the preprocessor variable
 *  #DEBUG# is set and if the given condition (in this case
 *  #n<dim# is violated.
 *
 *  If the exception was declared using the #DeclException0 (...)#
 *  macro, i.e. without any additional parameters, its name has
 *  nonetheless to be given with parentheses:
 *  #Assert (i>m, ExcSomewhat());#
 *
 *  \subsection{How it works internally}
 *  If the #DEBUG# preprocessor directive is set, the call #Assert (cond, exc);#
 *  is converted by the preprocessor into the sequence
 *  \begin{verbatim}
 *    if (!(cond))
 *      __IssueError_Assert (__FILE__,
 *                           __LINE__,
 *                           __PRETTY_FUNCTION__,
 *                           #cond,
 *                           #exc,
 *                           &exc);
 *  \end{verbatim}
 *  i.e. if the given condition is violated, then the file and
 *  line in which the exception occured as well as
 *  the condition itself and the call sequence of the
 *  exception object is transferred. Additionally an object
 *  of the form given by #exc# is created (this is normally an
 *  unnamed object like in #ExcDomain (n, dim)# of class
 *  #ExcDomain#) and transferred to the #__IssueError_Assert#
 *  function.
 *
 *  #__PRETTY__FUNCTION__# is a macro defined only by the GNU CC
 *  compiler and gives the name of the function. If another compiler
 *  is used, we set #__PRETTY_FUNCTION__ = "(unknown)"#.
 *
 *  In #__IssueError# the given data
 *  is transferred into the #exc# object by calling the
 *  #SetFields# function; after that, the general error info
 *  is printed onto #cerr# using the #PrintError# function of
 *  #exc# and finally the exception specific data is printed
 *  using the user defined function #PrintError# (which is
 *  normally created using the #DeclException (...)# macro
 *  family.
 *
 *  After printing all this information, #abort()# is called.
 *  This terminates the program, which is the right thing to do for
 *  this kind of error checking since it is used to detect programming
 *  errors rather than run-time errors; a program can, by definition,
 *  not recover from programming errors.
 *
 *  If the preprocessor variable #DEBUG# is not set, then nothing
 *  happens, i.e. the #Assert# macro is expanded to #{}#.
 *
 *  Sometimes, there is no useful condition for an exception other
 *  than that the program flow should not have reached a certain point,
 *  e.g. a #default# section of a #switch# statement. In this case,
 *  raise the exception by the following construct:
 *  \begin{verbatim}
 *    Assert (false, ExcInternalError());
 *  \end{verbatim}
 *
 *
 *  \section{Use of run-time exceptions}
 *
 *  For this mode, the standard #C++# #throw# and #catch# concept exists. We
 *  want to keep to this, but want to extend it a bit. In general, the
 *  structure is the same, i.e. you normally raise and exception by
 *  \begin{verbatim}
 *    if (!cond)
 *      throw ExcSomething();
 *  \end{verbatim}
 *  and catch it using the statement
 *  \begin{verbatim}
 *    try {
 *      do_something ();
 *    }
 *    catch (exception &e) {
 *      cerr << "Exception occured:" << endl
 *           << e.what ()
 *           << endl;
 *      do_something_to_reciver ();
 *    };
 *  \end{verbatim}
 *  #exception# is a standard #C++# class providing basic functionality for
 *  exceptions, such as the virtual function #what()# which returns some
 *  information on the exception itself. This information is useful if an
 *  exception can't be handled properly, in which case as precise a description
 *  as possible should be printed.
 *
 *  The problem here is that to get significant and useful information out
 *  of #what()#, it is necessary to overload this function in out exception
 *  class and call the #throw# operator with additional arguments to the
 *  exception class. The first thing, overloading the #what# function is
 *  done using the #DeclExceptionN# macros, but putting the right information,
 *  which is the same as explained above for the #Assert# expansion, requires
 *  some work if one would want to write it down each time:
 *  \begin{verbatim}
 *    if (!cond)
 *      {
 *        ExcSomething e(additional information);
 *        e.SetFields (__FILE__, __LINE__, __PRETTY_FUNCTION__,
 *                     "condition as a string",
 *                     "name of condition as a string");
 *        throw e;
 *      };
 *  \end{verbatim}
 *
 *  For this purpose, the macro #AssertThrow# was invented. It does mainly
 *  the same job as does the #Assert# macro, but it does not kill the
 *  program, it rather throws an exception as shown above. The mode of usage
 *  is
 *  \begin{verbatim}
 *    AssertThrow (cond, ExcSomething(additional information));
 *  \end{verbatim}
 *  The condition to be checked is incorporated into the macro in order to
 *  allow passing the violated condition as a string. The expansion of the
 *  #AssertThrow# macro is not affected by the #DEBUG# preprocessor variable.
 *
 *
 *  \section{Description of the #DeclExceptionN# macro family}
 *
 *  Declare an exception class without any additional parameters.
 *  There is a whole family of #DeclException?# macros
 *  where #?# is to be replaced by the number of additional
 *  parameters (0 to 5 presently).
 *  
 *  The syntax is as follows:
 *  \begin{verbatim}
 *    DeclException2 (ExcDomain,
 *                    int,
 *                    int,
 *                    << " i=" << arg1 << ", m=" << arg2);
 *  \end{verbatim}
 *  The first is the name of the exception class to be created.
 *  The next arguments are the types of the parameters (in this
 *  case there are two type names needed) and finally the output
 *  sequence with which you can print additional information.
 *  
 *  The syntax of the output sequence is a bit weird but gets
 *  clear once you see how this macro is defined:
 *  \begin{verbatim}
 *  class name : public ExceptionBase {
 *    public:
 *      name (const type1 a1, const type2 a2) :
 *                     arg1 (a1), arg2(a2) {};
 *      virtual void PrintInfo (ostream &out) const {
 *        out outsequence << endl;
 *      };
 *    private:
 *      type1 arg1;
 *      type2 arg2;
 *  };
 *  \end{verbatim}
 *   
 *  If declared as specified, you can later use this exception class
 *  in the following manner:
 *  \begin{verbatim}
 *    int i=5;
 *    int m=3;
 *    Assert (i<m, MyExc2(i,m));
 *  \end{verbatim}
 *  and the output if the condition fails will be
 *  \begin{verbatim}
 *    --------------------------------------------------------
 *    An error occurred in line <301> of file <exc-test.cc>.
 *    The violated condition was: 
 *      i<m
 *    The name and call sequence of the exception was:
 *      MyExc2(i,m)
 *    Additional Information: 
 *      i=5, m=3
 *    --------------------------------------------------------
 *  \end{verbatim}
 *  
 *  Obviously for the #DeclException0(name)# macro, no types and
 *  also no output sequence is allowed.
 *
 *
 *  @author Wolfgang Bangerth, November 1997, extensions 1998
 */
class ExceptionBase : public exception {
  public:
				     /**
				      * Default constructor.
				      */
    ExceptionBase ();
    
				     /**
				      *  The constructor takes the file in which the
				      *  error happened, the line and the violated
				      *  condition as well as the name of the
				      *  exception class as a #char*# as arguments.
				      */
    ExceptionBase (const char* f, const int l, const char *func,
		   const char* c, const char *e);

				     /**
				      *  Set the file name and line of where the
				      *  exception appeared as well as the violated
				      *  condition and the name of the exception as
				      *  a char pointer.
				      */
    void SetFields (const char *f,
		    const int   l,
		    const char *func,
		    const char *c,
		    const char *e);
    
				     /**
				      *  Print out the general part of the error
				      *  information.
				      */
    void PrintExcData (ostream &out) const;

				     /**
				      *  Print more specific information about the
				      *  exception which occured. Overload this
				      *  function in your own exception classes.
				      */
    virtual void PrintInfo (ostream &out) const;


/**
				      *  Function derived from the base class
				      *  which allows to pass information like
				      *  the line and name of the file where the
				      *  exception occured as well as user
				      *  information.
				      *
				      *  This function is mainly used when using
				      *  exceptions declared by the
				      *  #DeclException*# macros with the #throw#
				      *  mechanism or the #AssertThrow# macro.
				      */
    virtual const char * what () const;

  protected:
				     /**
				      * Name of the file this exception happen in.
				      */
    const char  *file;

				     /**
				      * Line number in this file.
				      */
    unsigned int line;

				     /**
				      * Name of the function, pretty printed.
				      */
    const char  *function;

				     /**
				      * The violated condition, as a string.
				      */
    const char  *cond;

				     /**
				      * Name of the exception and call sequence.
				      */
    const char  *exc;
};


/**
 *  This routine does the main work for the
 *  exception generation mechanism used in the
 *  #Assert# macro.
 *
 *  @see ExceptionBase
 */
template <class exc>
void __IssueError_Assert (const char *file,
			  int         line,
			  const char *function,
			  const char *cond,
			  const char *exc_name,
			  exc         e) {
				   // Fill the fields of the exception object
  e.SetFields (file, line, function, cond, exc_name);
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


/**
 *  This routine does the main work for the
 *  exception generation mechanism used in the
 *  #AssertThrow# macro.
 *
 *  @see ExceptionBase
 */
template <class exc>
void __IssueError_Throw (const char *file,
			 int         line,
			 const char *function,
			 const char *cond,
			 const char *exc_name,
			 exc         e) {
				   // Fill the fields of the exception object
  e.SetFields (file, line, function, cond, exc_name);
  throw e;
};


#ifdef DEBUG  ////////////////////////////////////////

/**
 *  This is the main routine in the exception mechanism for debug mode
 *  error checking. See the
 *  #ExceptionBase# class for more information.
 *
 *  @memo Assert that a certain condition is fulfilled, otherwise
 *   issue an error and abort the program.
 *  @see ExceptionBase
 *  @author Wolfgang Bangerth, November 1997, extensions 1998
 */
#define Assert(cond, exc)                                         \
  {                                                               \
    if (!(cond))                                                  \
      __IssueError_Assert (__FILE__,                              \
			   __LINE__,                              \
			   __PRETTY_FUNCTION__, #cond, #exc, exc);\
  }


#else        ////////////////////////////////////////

#define Assert(cond, exc)                     \
  { }
#endif      ////////////////////////////////////////


/**
 *  This is the main routine in the exception mechanism for run-time mode
 *  error checking. See the
 *  #ExceptionBase# class for more information.
 *
 *  @memo Assert that a certain condition is fulfilled, otherwise
 *   throw an exception
 *  @see ExceptionBase
 *  @author Wolfgang Bangerth, November 1997, extensions 1998
 */
#define AssertThrow(cond, exc)                                   \
  {                                                              \
    if (!(cond))                                                 \
      __IssueError_Throw (__FILE__,                              \
			  __LINE__,                              \
			  __PRETTY_FUNCTION__, #cond, #exc, exc);\
  }


/**
 * See the #ExceptionBase# class for a detailed description.
 *
 * @see ExceptionBase
 * @author Wolfgang Bangerth, November 1997
 */
#define DeclException0(Exception0)  \
class Exception0 :  public ExceptionBase {};

/**
  *  Declare an exception class with
  *  one additional parameter.
  *
  *  @see ExceptionBase
  */
#define DeclException1(Exception1, type1, outsequence)                \
class Exception1 : public ExceptionBase {                             \
  public:                                                             \
      Exception1 (const type1 a1) : arg1 (a1) {};                     \
      virtual void PrintInfo (ostream &out) const {                   \
        out outsequence << endl;                                      \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
}

/**
 *  Declare an exception class with
 *  two additional parameters.
 *
 *  @see ExceptionBase
 */
#define DeclException2(Exception2, type1, type2, outsequence)         \
class Exception2 : public ExceptionBase {                             \
  public:                                                             \
      Exception2 (const type1 a1, const type2 a2) :                   \
	      arg1 (a1), arg2(a2) {};                                 \
      virtual void PrintInfo (ostream &out) const {                   \
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
 *  @see ExceptionBase
 */
#define DeclException3(Exception3, type1, type2, type3, outsequence)  \
class Exception3 : public ExceptionBase {                             \
  public:                                                             \
      Exception3 (const type1 a1, const type2 a2, const type3 a3) :   \
	      arg1 (a1), arg2(a2), arg3(a3) {};                       \
      virtual void PrintInfo (ostream &out) const {                   \
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
 *  @see ExceptionBase
 */
#define DeclException4(Exception4, type1, type2, type3, type4, outsequence) \
class Exception4 : public ExceptionBase {                             \
  public:                                                             \
      Exception4 (const type1 a1, const type2 a2,                     \
	    const type3 a3, const type4 a4) :                         \
	      arg1 (a1), arg2(a2), arg3(a3), arg4(a4) {};             \
      virtual void PrintInfo (ostream &out) const {                   \
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
 *  @see ExceptionBase
 */
#define DeclException5(Exception5, type1, type2, type3, type4, type5, outsequence) \
class Exception5 : public ExceptionBase {                             \
  public:                                                             \
      Exception5 (const type1 a1, const type2 a2, const type3 a3,     \
	    const type4 a4, const type5 a5) :                         \
	      arg1 (a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {};   \
      virtual void PrintInfo (ostream &out) const {                   \
        out outsequence << endl;                                      \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
      const type2 arg2;                                               \
      const type3 arg3;                                               \
      const type4 arg4;                                               \
      const type5 arg5;                                               \
}


namespace StandardExceptions 
{
  DeclException0 (ExcDivideByZero);
  DeclException0 (ExcNotImplemented);
  DeclException0 (ExcInternalError);
  DeclException0 (ExcPureFunctionCalled);

  DeclException2 (ExcDimensionMismatch, int, int,
		  << "Dimension " << arg1 << " not equal to " << arg2);
  DeclException3 (ExcIndexRange, int, int, int,
		  << "Index " << arg1 << " is not in ["
		  << arg2 << "," << arg3 << "[");
  
};


using namespace StandardExceptions;


#endif
