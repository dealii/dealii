//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__exceptions_h
#define __deal2__exceptions_h

/**
 * @file
 * Here, the deal.II exception handling is located.
 */ 

#include <base/config.h>

#include <exception>

// we only need output streams, but older compilers did not provide
// them in a separate include file
#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif


/**
 *  This class should be used as a base class for
 *  all exception classes. Do not use its methods
 *  and variables directly since the interface
 *  and mechanism may be subject to change. Rather
 *  create new exception classes using the
 *  <tt>DeclException</tt> macro family.
 *
 *
 *  @section ExceptionBase General General overview of the exception handling mechanism in deal.II
 *
 *  The error handling mechanism in <tt>deal.II</tt> is generally used in two ways.
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
 *  program, rather they throw exceptions in the <tt>C++</tt> manner, allowing the
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
 *  additionally to the <tt>C++</tt> standard's <tt>exception</tt> class. Such a class
 *  is declared by the following lines of code:
 *   @code
 *     DeclException2 (ExcDomain, int, int,
 *                     << "Index= " << arg1 << "Upper Bound= " << arg2);
 *  @endcode
 *  This declares an exception class named <tt>ExcDomain</tt>, which
 *  has two variables as additional information (named
 *  <tt>arg1</tt> and <tt>arg2</tt> by default) and which outputs the
 *  given sequence (which is appended to an <tt>std::ostream</tt>
 *  variable's name, thus the weird syntax). There are
 *  other <tt>DeclExceptionN</tt> macros for exception classes
 *  with more or no parameters. It is proposed to let start
 *  the name of all exceptions by <tt>Exc...</tt> and to declare them locally to
 *  the class it is to be used in. Declaring exceptions globally is possible
 *  but pollutes the global namespace, is less readable and thus unnecessary.
 *
 *  Since exception classes are declared the same way for both modes of error
 *  checking, it is possible to use an exception declared through the
 *  <tt>DeclExceptionN(...)</tt> macro family in both modes; there is no need to
 *  declare different classes for each of these.
 *
 *
 *  @section ExceptionBaseUse Use of the debug mode exceptions
 *
 *  To use the exception mechanism for debug mode error checking, write lines
 *  like the following in your source code:
 *  @code
 *    Assert (n<dim, ExcDomain(n,dim));
 *  @endcode
 *  which by macro expansion does the following:
 *  @code
 *    #ifdef DEBUG
 *        if (!(cond))
 *              issue error of class ExcDomain(n,dim)
 *    #else
 *        do nothing
 *    #endif
 *  @endcode
 *  i.e. it issues an error only if the preprocessor variable
 *  <tt>DEBUG</tt> is set and if the given condition (in this case
 *  <tt>n < dim</tt> is violated).
 *
 *  If the exception was declared using the <tt>DeclException0 (...)</tt>
 *  macro, i.e. without any additional parameters, its name has
 *  nonetheless to be given with parentheses:
 *  <tt>Assert (i>m, ExcSomewhat());</tt>
 *
 *  @subsection ExceptionBaseInternal How it works internally
 *
 *  If the <tt>DEBUG</tt> preprocessor directive is set, the call <tt>Assert
 *  (cond, exc);</tt> is basically converted by the preprocessor into the
 *  sequence (note that function names and exact calling sequences may
 *  change over time, but the general principle remains the same)
 *  @code
 *    if (!(cond))
 *      deal_II_exceptions::internals::issue_error_assert_1
 *             (__FILE__,
 *              __LINE__,
 *              __PRETTY_FUNCTION__,
 *              #cond,
 *              #exc,
 *              &exc);
 *  @endcode
 *  i.e. if the given condition is violated, then the file and
 *  line in which the exception occured as well as
 *  the condition itself and the call sequence of the
 *  exception object is transferred. Additionally an object
 *  of the form given by <tt>exc</tt> is created (this is normally an
 *  unnamed object like in <tt>ExcDomain (n, dim)</tt> of class
 *  <tt>ExcDomain</tt>) and transferred to the deal_II_exceptions::internals::issue_error_assert_1()
 *  function.
 *
 *  <tt>__PRETTY__FUNCTION__</tt> is a macro defined by some compilers and
 *  gives the name of the function. If another compiler is used, we
 *  try to set this function to something reasonable, if the compiler
 *  provides us with that, and <tt>"(not available)"</tt> otherwise.
 *
 *  In <tt>__IssueError</tt> the given data
 *  is transferred into the <tt>exc</tt> object by calling the
 *  SetFields() function; after that, the general error info
 *  is printed onto <tt>std::cerr</tt> using the PrintError() function of
 *  <tt>exc</tt> and finally the exception specific data is printed
 *  using the user defined function PrintError() (which is
 *  normally created using the <tt>DeclException (...)</tt> macro
 *  family.
 *
 *  After printing all this information,
 *  deal_II_exceptions::internals::abort() is called (with one
 *  exception, see the end of this section). This terminates the
 *  program, which is the right thing to do for this kind of error
 *  checking since it is used to detect programming errors rather than
 *  run-time errors; a program can, by definition, not recover from
 *  programming errors.
 *
 *  If the preprocessor variable <tt>DEBUG</tt> is not set, then nothing
 *  happens, i.e. the <tt>Assert</tt> macro is expanded to <tt>{}</tt>.
 *
 *  Sometimes, there is no useful condition for an exception other
 *  than that the program flow should not have reached a certain point,
 *  e.g. a <tt>default</tt> section of a <tt>switch</tt> statement. In this case,
 *  raise the exception by the following construct:
 *  @verbatim
 *    Assert (false, ExcInternalError());
 *  @endverbatim
 *
 *  As mentioned above, the program is terminated once a call to
 *  <tt>Assert</tt> fails. However, there is one case where we do not want
 *  to do this, namely when a C++ exception is active. The usual case
 *  where this happens is that someone throws an exception through the
 *  <tt>AssertThrow</tt> mechanism (see below) which, while the stack is
 *  unwound, leads to the destruction of other objects in stack frames
 *  above. If other objects refer to the objects being thus destroyed,
 *  some destructors raise an exception through <tt>Assert</tt>. If we
 *  would abort the program then, we would only ever see the message
 *  that an object is being destroyed which is still referenced from
 *  somewhere, but we would never see the original exception that
 *  triggered this. (You can see it in the debugger by putting a break
 *  point on the function <tt>__throw</tt>, but you cannot see it from the
 *  program itself.) In that case, we use a C++ standard library
 *  function to detect the presence of another active exception and do
 *  not terminate the program to allow that the thrown exception
 *  propagates to some place where its message can be displayed.
 *
 *  Since it is common that one failed assertion leads to a whole
 *  chain of others, we only ever print the very first message. If the
 *  program is then aborted, that is no problem. If it is not (since a
 *  C++ exception is active), only the first is displayed and a
 *  message about suppressed follow-up messages is shown.
 *
 *
 *  @section ExceptionBaseUse2 Use of run-time exceptions
 *
 *  For this mode, the standard <tt>C++</tt> <tt>throw</tt> and <tt>catch</tt> concept exists. We
 *  want to keep to this, but want to extend it a bit. In general, the
 *  structure is the same, i.e. you normally raise and exception by
 *  @code
 *    if (!(cond))
 *      throw ExcSomething();
 *  @endcode
 *  and catch it using the statement
 *  @code
 *    try {
 *      do_something ();
 *    }
 *    catch (exception &e) {
 *      std::cerr << "Exception occured:" << std::endl
 *           << e.what ()
 *           << std::endl;
 *      do_something_to_reciver ();
 *    };
 *  @endcode
 *  <tt>exception</tt> is a standard <tt>C++</tt> class providing basic functionality for
 *  exceptions, such as the virtual function <tt>what()</tt> which returns some
 *  information on the exception itself. This information is useful if an
 *  exception can't be handled properly, in which case as precise a description
 *  as possible should be printed.
 *
 *  The problem here is that to get significant and useful information out
 *  of <tt>what()</tt>, it is necessary to overload this function in out exception
 *  class and call the <tt>throw</tt> operator with additional arguments to the
 *  exception class. The first thing, overloading the <tt>what</tt> function is
 *  done using the <tt>DeclExceptionN</tt> macros, but putting the right information,
 *  which is the same as explained above for the <tt>Assert</tt> expansion, requires
 *  some work if one would want to write it down each time:
 *  @code
 *    if (!(cond))
 *      {
 *        ExcSomething e(additional information);
 *        e.SetFields (__FILE__, __LINE__, __PRETTY_FUNCTION__,
 *                     "condition as a string",
 *                     "name of condition as a string");
 *        throw e;
 *      };
 *  @endcode
 *
 *  For this purpose, the macro <tt>AssertThrow</tt> was invented. It does mainly
 *  the same job as does the <tt>Assert</tt> macro, but it does not kill the
 *  program, it rather throws an exception as shown above. The mode of usage
 *  is
 *  @code
 *    AssertThrow (cond, ExcSomething(additional information));
 *  @endcode
 *  The condition to be checked is incorporated into the macro in order to
 *  allow passing the violated condition as a string. The expansion of the
 *  <tt>AssertThrow</tt> macro is not affected by the <tt>DEBUG</tt> preprocessor variable.
 *
 *
 *  @section ExceptionNMacros Description of the DeclExceptionN macro family
 *
 *  Declare an exception class without any additional parameters.
 *  There is a whole family of <tt>DeclException?</tt> macros
 *  where <tt>?</tt> is to be replaced by the number of additional
 *  parameters (0 to 5 presently).
 *  
 *  The syntax is as follows:
 *  @code
 *    DeclException2 (ExcDomain,
 *                    int,
 *                    int,
 *                    << " i=" << arg1 << ", m=" << arg2);
 *  @endcode
 *  The first is the name of the exception class to be created.
 *  The next arguments are the types of the parameters (in this
 *  case there are two type names needed) and finally the output
 *  sequence with which you can print additional information.
 *  
 *  The syntax of the output sequence is a bit weird but gets
 *  clear once you see how this macro is defined:
 *  @code
 *  class name : public ExceptionBase {
 *    public:
 *      name (const type1 a1, const type2 a2) :
 *                     arg1 (a1), arg2(a2) {};
 *      virtual void PrintInfo (std::ostream &out) const {
 *        out outsequence << std::endl;
 *      };
 *    private:
 *      type1 arg1;
 *      type2 arg2;
 *  };
 *  @endcode
 *   
 *  If declared as specified, you can later use this exception class
 *  in the following manner:
 *  @code
 *    int i=5;
 *    int m=3;
 *    Assert (i<m, MyExc2(i,m));
 *  @endcode
 *  and the output if the condition fails will be
 *  @verbatim
 *    --------------------------------------------------------
 *    An error occurred in line <301> of file <exc-test.cc>.
 *    The violated condition was: 
 *      i<m
 *    The name and call sequence of the exception was:
 *      MyExc2(i,m)
 *    Additional Information: 
 *      i=5, m=3
 *    --------------------------------------------------------
 *  @endverbatim
 *  
 *  Obviously for the <tt>DeclException0(name)</tt> macro, no types and
 *  also no output sequence is allowed.
 *
 *
 *  @author Wolfgang Bangerth, November 1997, extensions 1998
 */
class ExceptionBase : public std::exception
{
  public:
				     /**
				      * Default constructor.
				      */
    ExceptionBase ();
    
				     /**
				      *  The constructor takes the file in which the
				      *  error happened, the line and the violated
				      *  condition as well as the name of the
				      *  exception class as a <tt>char*</tt> as arguments.
				      */
    ExceptionBase (const char* f, const int l, const char *func,
		   const char* c, const char *e);

				     /**
				      * Destructor. Empty, but needed
				      * for the sake of exception
				      * specification, since the base
				      * class has this exception
				      * specification and the
				      * automatically generated
				      * destructor would have a
				      * different one due to member
				      * objects.
				      */
    virtual ~ExceptionBase () throw();
    
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
    void PrintExcData (std::ostream &out) const;

				     /**
				      *  Print out the stacktrace (if available)
                                      *  at the time the exception is created
				      */
    void PrintStackTrace (std::ostream &out) const;

				     /**
				      *  Print more specific information about the
				      *  exception which occured. Overload this
				      *  function in your own exception classes.
				      */
    virtual void PrintInfo (std::ostream &out) const;


				     /**
				      *  Function derived from the base class
				      *  which allows to pass information like
				      *  the line and name of the file where the
				      *  exception occured as well as user
				      *  information.
				      *
				      *  This function is mainly used when using
				      *  exceptions declared by the
				      *  <tt>DeclException*</tt> macros with the <tt>throw</tt>
				      *  mechanism or the <tt>AssertThrow</tt> macro.
				      */
    virtual const char * what () const throw ();

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
 * In this namespace functions in connection with the Assert and
 * AssertThrow mechanism are declared.
 */
namespace deal_II_exceptions
{

				   /**
				    * Set a string that is printed
				    * upon output of the message
				    * indicating a triggered
				    * <tt>Assert</tt> statement. This
				    * string, which is printed in
				    * addition to the usual output may
				    * indicate information that is
				    * otherwise not readily available
				    * unless we are using a
				    * debugger. For example, with
				    * distributed programs on cluster
				    * computers, the output of all
				    * processes is redirected to the
				    * same console window. In this
				    * case, it is convenient to set as
				    * additional name the name of the
				    * host on which the program runs,
				    * so that one can see in which
				    * instance of the program the
				    * exception occured.
				    *
				    * The string pointed to by the
				    * argument is copied, so needs not
				    * be stored after the call to this
				    * function.
				    *
				    * Previously set additional output
				    * is replaced by the argument
				    * given to this function.
				    */
  void set_additional_assert_output (const char * const p);
  
  
/**
 * The functions in this namespace are in connection with the Assert
 * and AssertThrow mechanism but are solely for internal purposes and
 * are not for use outside the exception handling and throwing
 * mechanism.
 */
  namespace internals
  {

/**
 *  This routine does the main work for the
 *  exception generation mechanism used in the
 *  <tt>Assert</tt> macro.
 *
 *  @ref ExceptionBase
 */
    void issue_error_assert (const char *file,
			     int         line,
			     const char *function,
			     const char *cond,
			     const char *exc_name,
			     ExceptionBase &e);
  

/**
 *  This routine does the main work for the
 *  exception generation mechanism used in the
 *  <tt>AssertThrow</tt> macro.
 *
 *  @ref ExceptionBase
 */
    template <class exc>
    void issue_error_throw (const char *file,
			    int         line,
			    const char *function,
			    const char *cond,
			    const char *exc_name,
			    exc         e)
    {
				       // Fill the fields of the
				       // exception object
      e.SetFields (file, line, function, cond, exc_name);
      throw e;
    }
    

/**
 *  Relay exceptions from the <tt>Assert</tt> macro to the
 *  <tt>__IssueError_Assert</tt> function. Used to convert the last
 *  argument from arbitrary type to ExceptionBase which is not
 *  possible inside the <tt>Assert</tt> macro due to syntactical
 *  difficulties in connection with the way we use the macro and the
 *  declaration of the exception classes.
 *
 *  @ref ExceptionBase
 */
    template <class exc>
    inline
    void issue_error_assert_1 (const char *file,
			       int         line,
			       const char *function,
			       const char *cond,
			       const char *exc_name,
			       exc         e)
    {
      issue_error_assert (file,line,function,cond,exc_name,e);
    }
    


/**
 * Abort the program. This function is used so that we need not
 * include <tt>cstdlib</tt> into this file since it is included into all
 * other files of the library and we would like to keep its include
 * list as short as possible.
 */
    void abort ();

  }
  
}



#ifdef DEBUG  ////////////////////////////////////////

/**
 * This is the main routine in the exception mechanism for debug mode
 * error checking. It asserts that a certain condition is fulfilled,
 * otherwise issues an error and aborts the program.
 *
 * See the ExceptionBase class for more information.
 *
 * @author Wolfgang Bangerth, November 1997, extensions 1998
 */
#define Assert(cond, exc)                                           \
  {                                                                 \
    if (!(cond))                                                    \
      deal_II_exceptions::internals::                               \
      issue_error_assert_1 (__FILE__,                               \
			     __LINE__,                              \
			     __PRETTY_FUNCTION__, #cond, #exc, exc);\
  }


#else        ////////////////////////////////////////

#define Assert(cond, exc)                     \
  { }
#endif      ////////////////////////////////////////



/**
 * This is the main routine in the exception mechanism for run-time
 * mode error checking. It assert that a certain condition is
 * fulfilled, otherwise issues an error and aborts the program.
 *
 * On some systems (we only know of DEC Alpha systems running under
 * OSF1 or Linux), the compiler fails to compile the <tt>AssertThrow</tt>
 * macro properly, yielding an internal compiler error. We detect this
 * at configure time. For these cases, the <tt>AssertThrow</tt> macro aborts
 * the program if the assertion is not satisfied. This, however,
 * happens in debug and optimized mode likewise.  Note that in these
 * cases, the meaning of a program changes. In particular, one cannot
 * catch exceptions thrown by <tt>AssertThrow</tt>, but we did not find
 * another way to work around this compiler bug.
 *
 * See the <tt>ExceptionBase</tt> class for more information.
 *
 * @ref ExceptionBase
 * @author Wolfgang Bangerth, November 1997, extensions 1998
 */
#ifndef DISABLE_ASSERT_THROW
#  ifndef HAVE_BUILTIN_EXPECT
#    define AssertThrow(cond, exc)                                   \
      {                                                              \
        if (!(cond))                                                 \
          deal_II_exceptions::internals::                            \
          issue_error_throw (__FILE__,                               \
		  	     __LINE__,                               \
			     __PRETTY_FUNCTION__, #cond, #exc, exc); \
      }
#  else // HAVE_BUILTIN_EXPECT
#    define AssertThrow(cond, exc)                                   \
      {                                                              \
        if (__builtin_expect(!(cond), false))                        \
          deal_II_exceptions::internals::                            \
          issue_error_throw (__FILE__,                               \
		  	     __LINE__,                               \
			     __PRETTY_FUNCTION__, #cond, #exc, exc); \
      }
#  endif
#else
#  define AssertThrow(cond, exc)                                    \
    {                                                               \
      if (!(cond))                                                  \
        deal_II_exceptions::internals::abort ();                    \
    }
#endif




/**
 * Declare an exception class derived from ExceptionBase without parameters.
 * @author Wolfgang Bangerth, November 1997
 */
#define DeclException0(Exception0)  \
class Exception0 :  public ExceptionBase {}



/**
  *  Declare an exception class derived from ExceptionBase with
  *  one additional parameter.
  */
#define DeclException1(Exception1, type1, outsequence)                \
class Exception1 : public ExceptionBase {                             \
  public:                                                             \
      Exception1 (const type1 a1) : arg1 (a1) {};                     \
      virtual ~Exception1 () throw () {};                             \
      virtual void PrintInfo (std::ostream &out) const {              \
        out outsequence << std::endl;                                 \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
}



/**
 *  Declare an exception class derived from ExceptionBase with
 *  two additional parameters.
 */
#define DeclException2(Exception2, type1, type2, outsequence)         \
class Exception2 : public ExceptionBase {                             \
  public:                                                             \
      Exception2 (const type1 a1, const type2 a2) :          \
	      arg1 (a1), arg2(a2) {};                                 \
      virtual ~Exception2 () throw () {};                             \
      virtual void PrintInfo (std::ostream &out) const {              \
        out outsequence << std::endl;                                 \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
      const type2 arg2;                                               \
}



/**
 *  Declare an exception class derived from ExceptionBase with
 *  three additional parameters.
 */
#define DeclException3(Exception3, type1, type2, type3, outsequence)  \
class Exception3 : public ExceptionBase {                             \
  public:                                                             \
      Exception3 (const type1 a1, const type2 a2, const type3 a3) : \
	      arg1 (a1), arg2(a2), arg3(a3) {};                       \
      virtual ~Exception3 () throw () {};                             \
      virtual void PrintInfo (std::ostream &out) const {              \
        out outsequence << std::endl;                                 \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
      const type2 arg2;                                               \
      const type3 arg3;                                               \
}



/**
 *  Declare an exception class derived from ExceptionBase with
 *  four additional parameters.
 */
#define DeclException4(Exception4, type1, type2, type3, type4, outsequence) \
class Exception4 : public ExceptionBase {                             \
  public:                                                             \
      Exception4 (const type1 a1, const type2 a2,                     \
	    const type3 a3, const type4 a4) :                \
	      arg1 (a1), arg2(a2), arg3(a3), arg4(a4) {};             \
      virtual ~Exception4 () throw () {};                             \
      virtual void PrintInfo (std::ostream &out) const {              \
        out outsequence << std::endl;                                 \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
      const type2 arg2;                                               \
      const type3 arg3;                                               \
      const type4 arg4;                                               \
}



/**
 *  Declare an exception class derived from ExceptionBase with
 *  five additional parameters.
 */
#define DeclException5(Exception5, type1, type2, type3, type4, type5, outsequence) \
class Exception5 : public ExceptionBase {                             \
  public:                                                             \
      Exception5 (const type1 a1, const type2 a2, const type3 a3,     \
	    const type4 a4, const type5 a5) :                \
	      arg1 (a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {};   \
      virtual ~Exception5 () throw () {};                             \
      virtual void PrintInfo (std::ostream &out) const {              \
        out outsequence << std::endl;                                 \
      };                                                              \
  private:                                                            \
      const type1 arg1;                                               \
      const type2 arg2;                                               \
      const type3 arg3;                                               \
      const type4 arg4;                                               \
      const type5 arg5;                                               \
}



/**
 * Declare some exceptions that occur over and over. This way, you can
 * simply use these exceptions, instead of having to declare them
 * locally in your class. The namespace in which these exceptions are
 * declared is later included into the global namespace by
 * @code
 * using namespace StandardExceptions;
 * @endcode
 */
namespace StandardExceptions 
{
				   /**
				    * @addtogroup Exceptions
				    */
				   //@{
  
				   /**
				    * Exception denoting a division by
				    * zero.
				    *
				    * @note Unfortunately, automatic
				    * detection of division by zero is
				    * very hardware dependent and
				    * requires severe hacking on some
				    * architectures. Therefore, this
				    * exception is only raised if the
				    * test s performed explicitly.
				    */
  DeclException0 (ExcDivideByZero);

				   /**
				    * Trying to allocate a new
				    * object failed due to lack of
				    * free memory.
				    */
  DeclException0 (ExcOutOfMemory);

				   /**
				    * An error occured reading or
				    * writing a file.
				    */
  DeclException0 (ExcIO);

				   /**
				    * An error occured opening the named file.
				    *
				    * The constructor takes a single
				    * argument of type <tt>char*</tt>
				    * naming the file.
				    */
  DeclException1 (ExcFileNotOpen,
		  char*,
		  << "Could not open file " << arg1);

				   /**
				    * Exception denoting a part of the
				    * library or application program
				    * that has not yet been
				    * implemented. In many cases, this
				    * only indicates that there wasn't
				    * much need for something yet, not
				    * that this is difficult to
				    * implement. It is therefore quite
				    * worth the effort to take a look
				    * at the corresponding place and
				    * see whether it can be
				    * implemented without too much
				    * effort.
				    */
  DeclException0 (ExcNotImplemented);

				   /**
				    * This exception usually indicates
				    * that some condition which the
				    * programmer thinks must be
				    * satisfied at a certain point in
				    * an algorithm, is not
				    * fulfilled. This might be due to
				    * some programming error above,
				    * due to changes to the algorithm
				    * that did not preserve this
				    * assertion, or due to assumptions
				    * the programmer made that are not
				    * valid at all (i.e. the exception
				    * is thrown although there is no
				    * error here). Within the library,
				    * this exception is most often
				    * used when we write some kind of
				    * complicated algorithm and are
				    * not yet sure whether we got it
				    * right; we then put in assertions
				    * after each part of the algorithm
				    * that check for some conditions
				    * that should hold there, and
				    * throw an exception if they do
				    * not.
				    *
				    * We usually leave in these
				    * assertions even after we are
				    * confident that the
				    * implementation is correct, since
				    * if someone later changes or
				    * extends the algorithm, these
				    * exceptions will indicate to him
				    * if he violates assumptions that
				    * are used later in the
				    * algorithm. Furthermore, it
				    * sometimes happens that an
				    * algorithm does not work in very
				    * rare corner cases. These cases
				    * will then be trapped sooner or
				    * later by the exception, so that
				    * the algorithm can then be fixed
				    * for these cases as well.
				    */
  DeclException0 (ExcInternalError);

				   /**
				    * This exception is used in
				    * functions that may not be called
				    * (i.e. in pure functions) but
				    * could not be declared pure since
				    * the class is intended to be used
				    * anyway, even though the
				    * respective function may only be
				    * called if a derived class is
				    * used.
				    */
  DeclException0 (ExcPureFunctionCalled);

                                   /**
				    * Used for constructors that are
				    * disabled.  Examples are copy
				    * constructors and assignment
				    * operators of large objects,
				    * which are only allowed for empty
				    * objects.
				    */
  DeclException0 (ExcInvalidConstructorCall);

				   /**
				    * This exception is used if some
				    * object is found uninitialized.
				    */
  DeclException0 (ExcNotInitialized);

				     /**
				      * The object is in a state not
				      * suitable for this operation.
				      */
    DeclException0 (ExcInvalidState);
    
				   /**
				    * This exception is raised if a
				    * functionality is not possible in
				    * the given dimension. Mostly used
				    * to throw function calls in 1d.
				    *
				    * The constructor takes a single
				    * <tt>int</tt>, denoting the
				    * dimension.
				    */
  DeclException1 (ExcImpossibleInDim,
		  int,
		  << "Impossible in " << arg1 << "d.");

				   /**
				    * A number is zero, but it should
				    * not be here.
				    */
  DeclException0(ExcZero);
  
				   /**
				    * The object should have been
				    * filled with something before
				    * this member function is called.
				    */
  DeclException0(ExcEmptyObject);
  
				   /**
				    * This exception is raised
				    * whenever the sizes of two
				    * objects were assumed to be
				    * equal, but were not.
				    *
				    * Parameters to the constructor
				    * are the first and second size,
				    * both of type <tt>int</tt>.
				    */
  DeclException2 (ExcDimensionMismatch,
		  int, int,
		  << "Dimension " << arg1 << " not equal to " << arg2);

				     /**
				      * The first dimension should be
				      * either equal to the second or
				      * the third, but it is neither.
				      */
    DeclException3 (ExcDimensionMismatch2,
		    int, int, int,
		    << "Dimension " << arg1 << " neither equal to " << arg2 << " nor to " << arg3);

				   /**
				    * This exception is one of the
				    * most often used ones, and
				    * indicates that an index is not
				    * within the expected range. For
				    * example, you might try to access
				    * an element of a vector which
				    * does not exist.
				    *
				    * The constructor takes three
				    * <tt>int</tt>, namely
				    * <ol>
				    * <li> the violating index
				    * <li> the lower bound
				    * <li> the upper bound plus one
				    * </ol>
				    */
  DeclException3 (ExcIndexRange, int, int, int,
		  << "Index " << arg1 << " is not in ["
		  << arg2 << "," << arg3 << "[");

				   /**
				    * A number is too small.
				    */
  DeclException2 (ExcLowerRange, int, int,
		  << "Number " << arg1
		  << " must be larger or equal " << arg2);
  
				   /**
				    * This exception indicates that
				    * the first argument should be an
				    * integer multiple of the second,
				    * but is not.
				    */
  DeclException2 (ExcNotMultiple, int, int,
		  << "Division " << arg1
		  << " by " << arg2
		  << " has remainder different from zero");
		  
				   /**
				    * This exception is thrown if the
				    * iterator you incremented or
				    * decremented was already at its
				    * final state.
				    */
  DeclException0 (ExcIteratorPastEnd);
  
				   /**
				    * This exception works around a
				    * design flaw in the
				    * <tt>DeclException0</tt> macro: that
				    * does not allow one to specify a
				    * message that is displayed when
				    * the exception is raised, as
				    * opposed to the other exceptions
				    * which allow to show a text along
				    * with the given parameters.
				    *
				    * When throwing this exception,
				    * you can give a message as a
				    * <tt>char*</tt> as argument to the
				    * exception that is then
				    * displayed.
				    */
  DeclException1 (ExcMessage, char*,
		  << arg1);

				   /**
				    * Exception used when running into
				    * functions that are only supported
				    * in a backward compatibility mode.
				    */
  DeclException1 (ExcCompatibility,
		  char*,
		  << "You are using a backward compatibility feature\n"
		  << "that you have disabled during configuration of\n"
		  << "the library by the --disable-compat="
		  << arg1 << " switch. You should either use an\n"
		  << "alternative function, or configure again without\n"
		  << "this switch and recompile the library.");

				   /**
				    * Some of our numerical classes
				    * allow for setting alll entries
				    * to zero using the assignment
				    * operator <tt>=</tt>.
				    *
				    * In many cases, this assignment
				    * operator makes sense <b>only</b>
				    * for the argument zero. In other
				    * cases, this exception is thrown.
				    */
  DeclException0 (ExcScalarAssignmentOnlyForZeroValue);

				   /**
				    * This function requires the BLAS
				    * library. Please reconfigure
				    * using the option
				    * <tt>--with-blas</tt> and check
				    * if it is actually included.
				    */
  DeclException0 (ExcNeedsBLAS);
  
				   /**
				    * This function requires the LAPACK
				    * library. Please reconfigure
				    * using the option
				    * <tt>--with-lapack</tt> and check
				    * if it is actually included.
				    */
  DeclException0 (ExcNeedsLAPACK);
  
				   /**
				    * This function requires the UMFPack
				    * library. Please reconfigure
				    * using the option
				    * <tt>--with-umfpack</tt> and check
				    * if it is actually included.
				    */
  DeclException0 (ExcNeedsUMFPACK);
  
				   /**
				    * This function requires the METIS
				    * library. Please reconfigure
				    * using the option
				    * <tt>--with-metis</tt> and check
				    * if it is actually included.
				    */
  DeclException0 (ExcNeedsMETIS);
  
				   /**
				    * This function requires the Petsc
				    * library. Please reconfigure
				    * using the option
				    * <tt>--with-petsc</tt> and check
				    * if it is actually included.
				    */
  DeclException0 (ExcNeedsPETSC);
				   /**
				    * A configuration option disabled
				    * this feature. In order to use
				    * it, you must reconfigure and
				    * recompile the libraries.
				    */
  DeclException1 (ExcDisabled, char*,
		  << "This feature was disabled by the "
		  "configuration option --disable-"
		  << arg1 << ". Reconfigure to use it!");
  
//@}
}

/*
 * Unfortunately, the following must be repeated for each library,
 * since we cannot have ifdefs in macros.
 */

/**
 * Assert support for the BLAS library
 */
#ifdef HAVE_LIBBLAS
#  define AssertBLAS {}
#else
#  define AssertBLAS Assert(false, ExcNeedsBLAS())
#endif


/**
 * Assert support for the LAPACK library
 */
#ifdef HAVE_LIBLAPACK
#  define AssertLAPACK {}
#else
#  define AssertLAPACK Assert(false, ExcNeedsLAPACK())
#endif


/**
 * Assert support for the UMFPACK library
 */
#ifdef HAVE_LIBUMFPACK
#  define AssertUMFPACK {}
#else
#  define AssertUMFPACK Assert(false, ExcNeedsUMFPACK())
#endif



using namespace StandardExceptions;


#endif
