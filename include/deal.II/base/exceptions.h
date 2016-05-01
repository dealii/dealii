// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#ifndef dealii__exceptions_h
#define dealii__exceptions_h

#include <deal.II/base/config.h>

#include <exception>
#include <string>
#include <ostream>

DEAL_II_NAMESPACE_OPEN


/**
 * This class is the base class for all exception classes. Do not use its
 * methods and variables directly since the interface and mechanism may be
 * subject to change. Rather create new exception classes using the
 * <tt>DeclException</tt> macro family.
 *
 * See the
 * @ref Exceptions
 * module for more details on this class and what can be done with classes
 * derived from it.
 *
 * @ingroup Exceptions
 * @author Wolfgang Bangerth, 1997, 1998, Matthias Maier, 2013
 */
class ExceptionBase : public std::exception
{
public:
  /**
   * Default constructor.
   */
  ExceptionBase ();

  /**
   * Copy constructor.
   */
  ExceptionBase (const ExceptionBase &exc);

  /**
   * Destructor.
   */
  virtual ~ExceptionBase () throw();

  /**
   * Set the file name and line of where the exception appeared as well as the
   * violated condition and the name of the exception as a char pointer. This
   * function also populates the stacktrace.
   */
  void set_fields (const char *file,
                   const int   line,
                   const char *function,
                   const char *cond,
                   const char *exc_name);


  /**
   * Override the standard function that returns the description of the error.
   */
  virtual const char *what() const throw();

  /**
   * Get exception name.
   */
  const char *get_exc_name() const;

  /**
   * Print out the general part of the error information.
   */
  void print_exc_data (std::ostream &out) const;

  /**
   * Print more specific information about the exception which occurred.
   * Overload this function in your own exception classes.
   */
  virtual void print_info (std::ostream &out) const;

  /**
   * Print a stacktrace, if one has been recorded previously, to the given
   * stream.
   */
  void print_stack_trace (std::ostream &out) const;

protected:
  /**
   * Name of the file this exception happens in.
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

  /**
   * A backtrace to the position where the problem happened, if the system
   * supports this.
   */
  mutable char **stacktrace;

  /**
   * The number of stacktrace frames that are stored in the previous variable.
   * Zero if the system does not support stack traces.
   */
  int n_stacktrace_frames;

#ifdef DEAL_II_HAVE_GLIBC_STACKTRACE
  /**
   * array of pointers that contains the raw stack trace
   */
  void *raw_stacktrace[25];
#endif

private:
  /**
   * Internal function that generates the c_string. Called by what().
   */
  void generate_message() const;

  /**
   * A pointer to the c_string that will be printed by what(). It is populated
   * by generate_message()
   */
  mutable std::string what_str;
};



/**
 * In this namespace, functions in connection with the Assert and AssertThrow
 * mechanism are declared.
 *
 * @ingroup Exceptions
 */
namespace deal_II_exceptions
{

  /**
   * Set a string that is printed upon output of the message indicating a
   * triggered <tt>Assert</tt> statement. This string, which is printed in
   * addition to the usual output may indicate information that is otherwise
   * not readily available unless we are using a debugger. For example, with
   * distributed programs on cluster computers, the output of all processes is
   * redirected to the same console window. In this case, it is convenient to
   * set as additional name the name of the host on which the program runs, so
   * that one can see in which instance of the program the exception occurred.
   *
   * The string pointed to by the argument is copied, so doesn't need to be
   * stored after the call to this function.
   *
   * Previously set additional output is replaced by the argument given to
   * this function.
   */
  void set_additional_assert_output (const char *const p);

  /**
   * Calling this function disables printing a stacktrace along with the other
   * output printed when an exception occurs. Most of the time, you will want
   * to see such a stacktrace; suppressing it, however, is useful if one wants
   * to compare the output of a program across different machines and systems,
   * since the stacktrace shows memory addresses and library names/paths that
   * depend on the exact setup of a machine.
   */
  void suppress_stacktrace_in_exceptions ();

  /**
   * Calling this function switches off the use of <tt>std::abort()</tt> when
   * an exception is created using the Assert() macro. Instead, the Exception
   * will be thrown using 'throw', so it can be caught if desired. Generally,
   * you want to abort the execution of a program when Assert() is called, but
   * it needs to be switched off if you want to log all exceptions created, or
   * if you want to test if an assertion is working correctly. This is done
   * for example in regression tests. Please note that some fatal errors will
   * still call abort(), e.g. when an exception is caught during exception
   * handling.
   */
  void disable_abort_on_exception ();

  /**
   * The functions in this namespace are in connection with the Assert and
   * AssertThrow mechanism but are solely for internal purposes and are not
   * for use outside the exception handling and throwing mechanism.
   *
   * @ingroup Exceptions
   */
  namespace internals
  {

    /**
     * Conditionally abort the program.
     *
     * Depending on whether disable_abort_on_exception was called, this
     * function either aborts the program flow by printing the error message
     * provided by @p exc and calling <tt>std::abort()</tt>, or throws @p exc
     * instead (if @p nothrow is set to <tt>false</tt>).
     *
     * If the boolean @p nothrow is set to true and disable_abort_on_exception
     * was called, the exception type is just printed to deallog and program
     * flow continues. This is useful if throwing an exception is prohibited
     * (e.g. in a destructor with <tt>noexcept(true)</tt> or
     * <tt>throw()</tt>).
     */
    void abort (const ExceptionBase &exc, bool nothrow = false);

    /**
     * An enum describing how to treat an exception in issue_error
     */
    enum ExceptionHandling
    {
      abort_on_exception,
      throw_on_exception,
      abort_nothrow_on_exception
    };

    /**
     * This routine does the main work for the exception generation mechanism
     * used in the <tt>Assert</tt> macro.
     *
     * @ref ExceptionBase
     */
    template <class exc>
    void issue_error (ExceptionHandling handling,
                      const char *file,
                      int         line,
                      const char *function,
                      const char *cond,
                      const char *exc_name,
                      exc         e)
    {
      // Fill the fields of the exception object
      e.set_fields (file, line, function, cond, exc_name);

      switch (handling)
        {
        case abort_on_exception:
          dealii::deal_II_exceptions::internals::abort(e);
          break;
        case abort_nothrow_on_exception:
          dealii::deal_II_exceptions::internals::abort(e, /*nothrow =*/ true);
          break;
        case throw_on_exception:
          throw e;
        }
    }

  } /*namespace internals*/

} /*namespace deal_II_exceptions*/



/**
 * This is the main routine in the exception mechanism for debug mode error
 * checking. It asserts that a certain condition is fulfilled, otherwise
 * issues an error and aborts the program.
 *
 * See the <tt>ExceptionBase</tt> class for more information.
 *
 * @note Active in DEBUG mode only
 * @ingroup Exceptions
 * @author Wolfgang Bangerth, 1997, 1998, Matthias Maier, 2013
 */
#ifdef DEBUG
#define Assert(cond, exc)                                                   \
  {                                                                           \
    if (!(cond))                                                              \
      ::dealii::deal_II_exceptions::internals::                               \
      issue_error(::dealii::deal_II_exceptions::internals::abort_on_exception,\
                  __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, #exc, exc); \
  }
#else
#define Assert(cond, exc)                                                   \
  {}
#endif



/**
 * A variant of the <tt>Assert</tt> macro above that exhibits the same runtime
 * behaviour as long as disable_abort_on_exception was not called.
 *
 * However, if disable_abort_on_exception was called, this macro merely prints
 * the exception that would be thrown to deallog and continues normally
 * without throwing an exception.
 *
 * See the <tt>ExceptionBase</tt> class for more information.
 *
 * @note Active in DEBUG mode only
 * @ingroup Exceptions
 * @author Wolfgang Bangerth, 1997, 1998, Matthias Maier, 2013
 */
#ifdef DEBUG
#define AssertNothrow(cond, exc)                                            \
  {                                                                           \
    if (!(cond))                                                              \
      ::dealii::deal_II_exceptions::internals::                               \
      issue_error(                                                            \
          ::dealii::deal_II_exceptions::internals::abort_nothrow_on_exception,  \
          __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, #exc, exc);           \
  }
#else
#define AssertNothrow(cond, exc)                                            \
  {}
#endif



/**
 * This is the main routine in the exception mechanism for run-time mode error
 * checking. It assert that a certain condition is fulfilled, otherwise issues
 * an error and aborts the program.
 *
 * See the <tt>ExceptionBase</tt> class for more information.
 *
 * @note Active in both DEBUG and RELEASE modes
 * @ingroup Exceptions
 * @author Wolfgang Bangerth, 1997, 1998, Matthias Maier, 2013
 */
#ifdef DEAL_II_HAVE_BUILTIN_EXPECT
#define AssertThrow(cond, exc)                                              \
  {                                                                           \
    if (__builtin_expect(!(cond), false))                                     \
      ::dealii::deal_II_exceptions::internals::                               \
      issue_error(::dealii::deal_II_exceptions::internals::throw_on_exception,\
                  __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, #exc, exc); \
  }
#else /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/
#define AssertThrow(cond, exc)                                              \
  {                                                                           \
    if (!(cond))                                                              \
      ::dealii::deal_II_exceptions::internals::                               \
      issue_error(::dealii::deal_II_exceptions::internals::throw_on_exception,\
                  __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, #exc, exc); \
  }
#endif /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/



#ifndef DOXYGEN

/**
 * Declare an exception class derived from ExceptionBase without parameters.
 *
 * @author Wolfgang Bangerth, November 1997
 * @ingroup Exceptions
 */
#define DeclException0(Exception0)                                        \
  class Exception0 :  public dealii::ExceptionBase {}


/**
 * Declare an exception class derived from ExceptionBase that can take one
 * runtime argument, but if none is given in the place where you want to throw
 * the exception, it simply reverts to the default text provided when
 * declaring the exception class through this macro.
 *
 * @ingroup Exceptions
 */
#define DeclExceptionMsg(Exception, defaulttext)                          \
  class Exception :  public dealii::ExceptionBase                         \
  {                                                                       \
  public:                                                                 \
    Exception (const std::string &msg = defaulttext) : arg (msg) {}       \
    virtual ~Exception () throw () {}                                     \
    virtual void print_info (std::ostream &out) const {                   \
      out << arg << std::endl;                                            \
    }                                                                     \
  private:                                                                \
    const std::string arg;                                                \
  }

/**
 * Declare an exception class derived from ExceptionBase with one additional
 * parameter.
 *
 * @ingroup Exceptions
 */
#define DeclException1(Exception1, type1, outsequence)                    \
  class Exception1 : public dealii::ExceptionBase {                       \
  public:                                                                 \
    Exception1 (const type1 a1) : arg1 (a1) {}                            \
    virtual ~Exception1 () throw () {}                                    \
    virtual void print_info (std::ostream &out) const {                   \
      out outsequence << std::endl;                                       \
    }                                                                     \
  private:                                                                \
    const type1 arg1;                                                     \
  }


/**
 * Declare an exception class derived from ExceptionBase with two additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DeclException2(Exception2, type1, type2, outsequence)             \
  class Exception2 : public dealii::ExceptionBase {                       \
  public:                                                                 \
    Exception2 (const type1 a1, const type2 a2) :                         \
      arg1 (a1), arg2(a2) {}                                              \
    virtual ~Exception2 () throw () {}                                    \
    virtual void print_info (std::ostream &out) const {                   \
      out outsequence << std::endl;                                       \
    }                                                                     \
  private:                                                                \
    const type1 arg1;                                                     \
    const type2 arg2;                                                     \
  }


/**
 * Declare an exception class derived from ExceptionBase with three additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DeclException3(Exception3, type1, type2, type3, outsequence)      \
  class Exception3 : public dealii::ExceptionBase {                       \
  public:                                                                 \
    Exception3 (const type1 a1, const type2 a2, const type3 a3) :         \
      arg1 (a1), arg2(a2), arg3(a3) {}                                    \
    virtual ~Exception3 () throw () {}                                    \
    virtual void print_info (std::ostream &out) const {                   \
      out outsequence << std::endl;                                       \
    }                                                                     \
  private:                                                                \
    const type1 arg1;                                                     \
    const type2 arg2;                                                     \
    const type3 arg3;                                                     \
  }


/**
 * Declare an exception class derived from ExceptionBase with four additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DeclException4(Exception4, type1, type2, type3, type4, outsequence) \
  class Exception4 : public dealii::ExceptionBase {                       \
  public:                                                                 \
    Exception4 (const type1 a1, const type2 a2,                           \
                const type3 a3, const type4 a4) :                         \
      arg1 (a1), arg2(a2), arg3(a3), arg4(a4) {}                          \
    virtual ~Exception4 () throw () {}                                    \
    virtual void print_info (std::ostream &out) const {                   \
      out outsequence << std::endl;                                       \
    }                                                                     \
  private:                                                                \
    const type1 arg1;                                                     \
    const type2 arg2;                                                     \
    const type3 arg3;                                                     \
    const type4 arg4;                                                     \
  }


/**
 * Declare an exception class derived from ExceptionBase with five additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DeclException5(Exception5, type1, type2, type3, type4, type5, outsequence) \
  class Exception5 : public dealii::ExceptionBase {                       \
  public:                                                                 \
    Exception5 (const type1 a1, const type2 a2, const type3 a3,           \
                const type4 a4, const type5 a5) :                         \
      arg1 (a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {}                \
    virtual ~Exception5 () throw () {}                                    \
    virtual void print_info (std::ostream &out) const {                   \
      out outsequence << std::endl;                                       \
    }                                                                     \
  private:                                                                \
    const type1 arg1;                                                     \
    const type2 arg2;                                                     \
    const type3 arg3;                                                     \
    const type4 arg4;                                                     \
    const type5 arg5;                                                     \
  }

#else /*ifndef DOXYGEN*/

// Dummy definitions for doxygen:

/**
 * Declare an exception class derived from ExceptionBase without parameters.
 *
 * @author Wolfgang Bangerth, November 1997
 * @ingroup Exceptions
 */
#define DeclException0(Exception0)                                        \
  /** @ingroup Exceptions */ \
  static dealii::ExceptionBase& Exception0 ()

/**
 * Declare an exception class derived from ExceptionBase that can take one
 * runtime argument, but if none is given in the place where you want to throw
 * the exception, it simply reverts to the default text provided when
 * declaring the exception class through this macro.
 *
 * @ingroup Exceptions
 */
#define DeclExceptionMsg(Exception, defaulttext)                          \
  /** @ingroup Exceptions */ \
  /** @note The error text is: @code defaulttext @endcode */  \
  static dealii::ExceptionBase& Exception ()

/**
 * Declare an exception class derived from ExceptionBase with one additional
 * parameter.
 *
 * @ingroup Exceptions
 */
#define DeclException1(Exception1, type1, outsequence)                    \
  /** @ingroup Exceptions */ \
  /** @note The error text is: @code outsequence @endcode */              \
  static dealii::ExceptionBase& Exception1 (type1 arg1)


/**
 * Declare an exception class derived from ExceptionBase with two additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DeclException2(Exception2, type1, type2, outsequence)             \
  /** @ingroup Exceptions */ \
  /** @note The error text is: @code outsequence @endcode */              \
  static dealii::ExceptionBase& Exception2 (type1 arg1, type2 arg2)


/**
 * Declare an exception class derived from ExceptionBase with three additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DeclException3(Exception3, type1, type2, type3, outsequence)      \
  /** @ingroup Exceptions */ \
  /** @note The error text is: @code outsequence @endcode */              \
  static dealii::ExceptionBase& Exception3 (type1 arg1, type2 arg2, type3 arg3)


/**
 * Declare an exception class derived from ExceptionBase with four additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DeclException4(Exception4, type1, type2, type3, type4, outsequence) \
  /** @ingroup Exceptions */ \
  /** @note The error text is: @code outsequence @endcode */                \
  static dealii::ExceptionBase& Exception4 (type1 arg1, type2 arg2, type3 arg3, type4 arg4)


/**
 * Declare an exception class derived from ExceptionBase with five additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DeclException5(Exception5, type1, type2, type3, type4, type5, outsequence) \
  /** @ingroup Exceptions */ \
  /** @note The error text is: @code outsequence @endcode */                       \
  static dealii::ExceptionBase& Exception5 (type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5)

#endif /*ifndef DOXYGEN*/


/**
 * Declare some exceptions that occur over and over. This way, you can simply
 * use these exceptions, instead of having to declare them locally in your
 * class. The namespace in which these exceptions are declared is later
 * included into the global namespace by
 * @code
 * using namespace StandardExceptions;
 * @endcode
 *
 * @ingroup Exceptions
 */
namespace StandardExceptions
{
  /**
   * @addtogroup Exceptions
   */
  //@{

  /**
   * Exception denoting a division by zero.
   */
  DeclExceptionMsg (ExcDivideByZero,
                    "A piece of code is attempting a division by zero. This is "
                    "likely going to lead to results that make no sense.");

  /**
   * Exception raised if a number is not finite.
   *
   * This exception should be used to catch infinite or not a number results
   * of arithmetic operations that do not result from a division by zero (use
   * ExcDivideByZero for those).
   *
   * The exception uses std::complex as its argument to ensure that we can use
   * it for all scalar arguments (real or complex-valued).
   */
  DeclException1 (ExcNumberNotFinite,
                  std::complex<double>,
                  << "In a significant number of places, deal.II checks that some intermediate "
                  << "value is a finite number (as opposed to plus or minus infinity, or "
                  << "NaN/Not a Number). In the current function, we encountered a number "
                  << "that is not finite (its value is " << arg1 << " and therefore "
                  << "violates the current assertion.\n\n"
                  << "This may be due to the fact that some operation in this function "
                  << "created such a value, or because one of the arguments you passed "
                  << "to the function already had this value from some previous "
                  << "operation. In the latter case, this function only triggered the "
                  << "error but may not actually be responsible for the computation of "
                  << "the number that is not finite.\n\n"
                  << "There are two common cases where this situation happens. First, your "
                  << "code (or something in deal.II) divides by zero in a place where this "
                  << "should not happen. Or, you are trying to solve a linear system "
                  << "with an unsuitable solver (such as an indefinite or non-symmetric "
                  << "linear system using a Conjugate Gradient solver); such attempts "
                  << "oftentimes yield an operation somewhere that tries to divide "
                  << "by zero or take the square root of a negative value.\n\n"
                  << "In any case, when trying to find the source of the error, "
                  << "recall that the location where you are getting this error is "
                  << "simply the first place in the program where there is a check "
                  << "that a number (e.g., an element of a solution vector) is in fact "
                  << "finite, but that the actual error that computed the number "
                  << "may have happened far earlier. To find this location, you "
                  << "may want to add checks for finiteness in places of your "
                  << "program visited before the place where this error is produced."
                  << "One way to check for finiteness is to use the 'AssertIsFinite' "
                  << "macro.");

  /**
   * Trying to allocate a new object failed due to lack of free memory.
   */
  DeclExceptionMsg (ExcOutOfMemory,
                    "Your program tried to allocate some memory but this "
                    "allocation failed. Typically, this either means that "
                    "you simply do not have enough memory in your system, "
                    "or that you are (erroneously) trying to allocate "
                    "a chunk of memory that is simply beyond all reasonable "
                    "size, for example because the size of the object has "
                    "been computed incorrectly.");

  /**
   * A memory handler reached a point where all allocated objects should have
   * been released. Since this exception is thrown, some were still allocated.
   */
  DeclException1 (ExcMemoryLeak, int,
                  << "Destroying memory handler while " << arg1
                  << " objects are still allocated");

  /**
   * An error occurred reading or writing a file.
   */
  DeclExceptionMsg (ExcIO,
                    "An input/output error has occurred. There are a number of "
                    "reasons why this may be happening, both for reading and "
                    "writing operations."
                    "\n\n"
                    "If this happens during an operation that tries to read "
                    "data: First, you may be "
                    "trying to read from a file that doesn't exist or that is "
                    "not readable given its file permissions. Second, deal.II "
                    "uses this error at times if it tries to "
                    "read information from a file but where the information "
                    "in the file does not correspond to the expected format. "
                    "An example would be a truncated file, or a mesh file "
                    "that contains not only sections that describe the "
                    "vertices and cells, but also sections for additional "
                    "data that deal.II does not understand."
                    "\n\n"
                    "If this happens during an operation that tries to write "
                    "data: you may be trying to write to a file to which file "
                    "or directory permissions do not allow you to write. A "
                    "typical example is where you specify an output file in "
                    "a directory that does not exist.");

  /**
   * An error occurred opening the named file.
   *
   * The constructor takes a single argument of type <tt>char*</tt> naming the
   * file.
   */
  DeclException1 (ExcFileNotOpen,
                  char *,
                  << "Could not open file " << arg1);

  /**
   * Exception denoting a part of the library or application program that has
   * not yet been implemented. In many cases, this only indicates that there
   * wasn't much need for something yet, not that this is difficult to
   * implement. It is therefore quite worth the effort to take a look at the
   * corresponding place and see whether it can be implemented without too
   * much effort.
   */
  DeclExceptionMsg (ExcNotImplemented,
                    "You are trying to use functionality in deal.II that is "
                    "currently not implemented. In many cases, this indicates "
                    "that there simply didn't appear much of a need for it, or "
                    "that the author of the original code did not have the "
                    "time to implement a particular case. If you hit this "
                    "exception, it is therefore worth the time to look into "
                    "the code to find out whether you may be able to "
                    "implement the missing functionality. If you do, please "
                    "consider providing a patch to the deal.II development "
                    "sources (see the deal.II website on how to contribute).");

  /**
   * This exception usually indicates that some condition which the programmer
   * thinks must be satisfied at a certain point in an algorithm, is not
   * fulfilled. This might be due to some programming error above, due to
   * changes to the algorithm that did not preserve this assertion, or due to
   * assumptions the programmer made that are not valid at all (i.e. the
   * exception is thrown although there is no error here). Within the library,
   * this exception is most often used when we write some kind of complicated
   * algorithm and are not yet sure whether we got it right; we then put in
   * assertions after each part of the algorithm that check for some
   * conditions that should hold there, and throw an exception if they do not.
   *
   * We usually leave in these assertions even after we are confident that the
   * implementation is correct, since if someone later changes or extends the
   * algorithm, these exceptions will indicate to him if he violates
   * assumptions that are used later in the algorithm. Furthermore, it
   * sometimes happens that an algorithm does not work in very rare corner
   * cases. These cases will then be trapped sooner or later by the exception,
   * so that the algorithm can then be fixed for these cases as well.
   */
  DeclExceptionMsg (ExcInternalError,
                    "This exception -- which is used in many places in the "
                    "library -- usually indicates that some condition which "
                    "the author of the code thought must be satisfied at a "
                    "certain point in an algorithm, is not fulfilled. An "
                    "example would be that the first part of an algorithm "
                    "sorts elements of an array in ascending order, and "
                    "a second part of the algorithm later encounters an "
                    "an element that is not larger than the previous one."
                    "\n\n"
                    "There is usually not very much you can do if you "
                    "encounter such an exception since it indicates an error "
                    "in deal.II, not in your own program. Try to come up with "
                    "the smallest possible program that still demonstrates "
                    "the error and contact the deal.II mailing lists with it "
                    "to obtain help.");

  /**
   * This exception is used in functions that may not be called (i.e. in pure
   * functions) but could not be declared pure since the class is intended to
   * be used anyway, even though the respective function may only be called if
   * a derived class is used.
   */
  DeclExceptionMsg (ExcPureFunctionCalled,
                    "You (or a place in the library) are trying to call a "
                    "function that is declared as a virtual function in a "
                    "base class but that has not been overridden in your "
                    "derived class."
                    "\n\n"
                    "This exception happens in cases where the base class "
                    "cannot provide a useful default implementation for "
                    "the virtual function, but where we also do not want "
                    "to mark the function as abstract (i.e., with '=0' at the end) "
                    "because the function is not essential to the class in many "
                    "contexts. In cases like this, the base class provides "
                    "a dummy implementation that makes the compiler happy, but "
                    "that then throws the current exception."
                    "\n\n"
                    "A concrete example would be the 'Function' class. It declares "
                    "the existence of 'value()' and 'gradient()' member functions, "
                    "and both are marked as 'virtual'. Derived classes have to "
                    "override these functions for the values and gradients of a "
                    "particular function. On the other hand, not every function "
                    "has a gradient, and even for those that do, not every program "
                    "actually needs to evaluate it. Consequently, there is no "
                    "*requirement* that a derived class actually override the "
                    "'gradient()' function (as there would be had it been marked "
                    "as abstract). But, since the base class cannot know how to "
                    "compute the gradient, if a derived class does not override "
                    "the 'gradient()' function and it is called anyway, then the "
                    "default implementation in the base class will simply throw "
                    "an exception."
                    "\n\n"
                    "The exception you see is what happens in cases such as the "
                    "one just illustrated. To fix the problem, you need to "
                    "investigate whether the function being called should indeed have "
                    "been called; if the answer is 'yes', then you need to "
                    "implement the missing override in your class.");

  /**
   * Used for constructors that are disabled. Examples are copy constructors
   * and assignment operators of large objects, which are only allowed for
   * empty objects.
   */
  DeclException0 (ExcInvalidConstructorCall);

  /**
   * This exception is used if some object is found uninitialized.
   */
  DeclException0 (ExcNotInitialized);

  /**
   * The object is in a state not suitable for this operation.
   */
  DeclException0 (ExcInvalidState);

  /**
   * This exception is raised if a functionality is not possible in the given
   * dimension. Mostly used to throw function calls in 1d.
   *
   * The constructor takes a single <tt>int</tt>, denoting the dimension.
   */
  DeclException1 (ExcImpossibleInDim,
                  int,
                  << "You are trying to execute functionality that is "
                  << "impossible in " << arg1
                  << "d or simply does not make any sense.");

  /**
   * This exception is raised if a functionality is not possible in the given
   * combination of dimension and space-dimension.
   *
   * The constructor takes two <tt>int</tt>, denoting the dimension and the
   * space dimension.
   */
  DeclException2 (ExcImpossibleInDimSpacedim,
                  int, int,
                  << "You are trying to execute functionality that is "
                  << "impossible in dimensions <" << arg1 << "," << arg2
                  << "> or simply does not make any sense.");


  /**
   * A number is zero, but it should not be here.
   */
  DeclExceptionMsg(ExcZero,
                   "In a check in the code, deal.II encountered a zero in "
                   "a place where this does not make sense. See the condition "
                   "that was being checked and that is printed further up "
                   "in the error message to get more information on what "
                   "the erroneous zero corresponds to.");

  /**
   * The object should have been filled with something before this member
   * function is called.
   */
  DeclExceptionMsg(ExcEmptyObject,
                   "The object you are trying to access is empty but it makes "
                   "no sense to attempt the operation you are trying on an "
                   "empty object.");

  /**
   * This exception is raised whenever the sizes of two objects were assumed
   * to be equal, but were not.
   *
   * Parameters to the constructor are the first and second size, both of type
   * <tt>int</tt>.
   */
  DeclException2 (ExcDimensionMismatch,
                  std::size_t, std::size_t,
                  << "Dimension " << arg1 << " not equal to " << arg2);

  /**
   * The first dimension should be either equal to the second or the third,
   * but it is neither.
   */
  DeclException3 (ExcDimensionMismatch2,
                  int, int, int,
                  << "Dimension " << arg1 << " neither equal to " << arg2
                  << " nor to " << arg3);

  /**
   * This exception indicates that an index is not within the expected range.
   * For example, it may be that you are trying to access an element of a
   * vector which does not exist.
   *
   * The constructor takes three <tt>int</tt> arguments, namely
   * <ol>
   * <li> the violating index
   * <li> the lower bound
   * <li> the upper bound plus one
   * </ol>
   */
  DeclException3 (ExcIndexRange,
                  int, int, int,
                  << "Index " << arg1 << " is not in the half-open range [" << arg2 << ","
                  << arg3 << ")."
                  << (arg2==arg3 ?
                      " In the current case, this half-open range is in fact empty, "
                      "suggesting that you are accessing an element of an empty "
                      "collection such as a vector that has not been set to the "
                      "correct size."
                      :
                      ""));

  /**
   * This exception indicates that an index is not within the expected range.
   * For example, it may be that you are trying to access an element of a
   * vector which does not exist.
   *
   * The constructor takes three <tt>int</tt> arguments, namely
   * <ol>
   * <li> the violating index
   * <li> the lower bound
   * <li> the upper bound plus one
   * </ol>
   *
   * This generic exception differs from ExcIndexRange by allowing to specify
   * the type of indices.
   */
  template <typename T>
  DeclException3 (ExcIndexRangeType,
                  T,T,T,
                  << "Index " << arg1 << " is not in the half-open range [" << arg2 << ","
                  << arg3 << ")."
                  << (arg2==arg3 ?
                      " In the current case, this half-open range is in fact empty, "
                      "suggesting that you are accessing an element of an empty "
                      "collection such as a vector that has not been set to the "
                      "correct size."
                      :
                      ""));

  /**
   * A number is too small.
   */
  DeclException2 (ExcLowerRange,
                  int, int,
                  << "Number " << arg1 << " must be larger than or equal "
                  << arg2 << ").");

  /**
   * A generic exception definition for the ExcLowerRange above.
   */
  template <typename T>
  DeclException2 (ExcLowerRangeType,
                  T, T,
                  << "Number " << arg1 << " must be larger than or equal "
                  << arg2 << ".");

  /**
   * This exception indicates that the first argument should be an integer
   * multiple of the second, but is not.
   */
  DeclException2 (ExcNotMultiple,
                  int, int,
                  << "Division " << arg1 << " by " << arg2
                  << " has remainder different from zero");

  /**
   * This exception is thrown if the iterator you access has corrupted data.
   * It might for instance be, that the container it refers does not have an
   * entry at the point the iterator refers.
   *
   * Typically, this will be an internal error of deal.II, because the
   * increment and decrement operators should never yield an invalid iterator.
   */
  DeclExceptionMsg (ExcInvalidIterator,
                    "You are trying to use an iterator, but the iterator is "
                    "in an invalid state. This may indicate that the iterator "
                    "object has not been initialized, or that it has been "
                    "moved beyond the end of the range of valid elements.");

  /**
   * This exception is thrown if the iterator you incremented or decremented
   * was already at its final state.
   */
  DeclExceptionMsg (ExcIteratorPastEnd,
                    "You are trying to use an iterator, but the iterator is "
                    "pointing past the end of the range of valid elements. "
                    "It is not valid to dereference the iterator in this "
                    "case.");

  /**
   * This exception works around a design flaw in the <tt>DeclException0</tt>
   * macro: exceptions declared through DeclException0 do not allow one to
   * specify a message that is displayed when the exception is raised, as
   * opposed to the other exceptions which allow to show a text along with the
   * given parameters.
   *
   * When throwing this exception, you can give a message as a
   * <tt>std::string</tt> as argument to the exception that is then displayed.
   * The argument can, of course, be constructed at run-time, for example
   * including the name of a file that can't be opened, or any other text you
   * may want to assemble from different pieces.
   */
  DeclException1 (ExcMessage,
                  std::string,
                  << arg1);

  /**
   * Parallel vectors with ghost elements are read-only vectors.
   */
  DeclExceptionMsg (ExcGhostsPresent,
                    "You are trying an operation on a vector that is only "
                    "allowed if the vector has no ghost elements, but the "
                    "vector you are operating on does have ghost elements. "
                    "Specifically, vectors with ghost elements are read-only "
                    "and cannot appear in operations that write into these "
                    "vectors."
                    "\n\n"
                    "See the glossary entry on 'Ghosted vectors' for more "
                    "information.");

  /**
   * Some of our numerical classes allow for setting all entries to zero using
   * the assignment operator <tt>=</tt>.
   *
   * In many cases, this assignment operator makes sense <b>only</b> for the
   * argument zero. In other cases, this exception is thrown.
   */
  DeclExceptionMsg (ExcScalarAssignmentOnlyForZeroValue,
                    "You are trying an operation of the form 'vector=s' with "
                    "a nonzero scalar value 's'. However, such assignments "
                    "are only allowed if the right hand side is zero.");

  /**
   * This function requires support for the LAPACK library.
   */
  DeclException0 (ExcNeedsLAPACK);

  /**
   * This function requires support for the NetCDF library.
   */
  DeclException0 (ExcNeedsNetCDF);

  /**
   * This function requires support for the FunctionParser library.
   */
  DeclException0 (ExcNeedsFunctionparser);


//@}
} /*namespace StandardExceptions*/


/**
 * Special assertion for dimension mismatch.
 *
 * Since this is used very often and always repeats the arguments, we
 * introduce this special assertion for ExcDimensionMismatch in order to keep
 * the user codes shorter.
 *
 * @ingroup Exceptions
 * @author Guido Kanschat 2007
 */
#define AssertDimension(dim1,dim2) Assert((dim1) == (dim2),       \
                                          dealii::ExcDimensionMismatch((dim1),(dim2)))


/**
 * An assertion that tests whether <tt>vec</tt> has size <tt>dim1</tt>, and
 * each entry of the vector is itself an array that has the size <tt>dim2</tt>.
 *
 * @ingroup Exceptions
 * @author Guido Kanschat 2010
 */
#define AssertVectorVectorDimension(vec,dim1,dim2) AssertDimension((vec).size(), (dim1)) \
  for (unsigned int i=0;i<dim1;++i) { AssertDimension((vec)[i].size(), (dim2)); }


/**
 * An assertion that tests that a given index is within the half-open
 * range <code>[0,range)</code>. It throws an exception object
 * <code>ExcIndexRange(index,0,range)</code> if the assertion
 * fails.
 *
 * @ingroup Exceptions
 * @author Guido Kanschat 2007
 */
#define AssertIndexRange(index,range) Assert((index) < (range), \
                                             dealii::ExcIndexRange((index),0,(range)))

/**
 * An assertion that tests that a given index is within the half-open
 * range <code>[0,range)</code>. It throws an exception object
 * <code>ExcIndexRange(index,0,range)</code> if the assertion
 * fails.
 *
 * This variation of the AssertIndexRange assertion is used for indices of type
 * types::global_dof_index for which we need to make special accommodations because
 * they may not fit into the regular 32-bit integer indices used in AssertIndexRange.
 *
 * @ingroup Exceptions
 */
#define AssertGlobalIndexRange(index,range) Assert((index) < (range), \
                                                   dealii::ExcIndexRange<types::global_dof_index>((index),0,(range)))

/**
 * An assertion that checks whether a number is finite or not. We explicitly
 * cast the number to std::complex to match the signature of the exception
 * (see there for an explanation of why we use std::complex at all) and to
 * satisfy the fact that std::complex has no implicit conversions.
 *
 * @ingroup Exceptions
 * @author Wolfgang Bangerth, 2015
 */
#define AssertIsFinite(number) Assert(dealii::numbers::is_finite(number), \
                                      dealii::ExcNumberNotFinite(std::complex<double>(number)))

using namespace StandardExceptions;

DEAL_II_NAMESPACE_CLOSE

#endif
