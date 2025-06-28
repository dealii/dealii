// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_exceptions_h
#define dealii_exceptions_h

#include <deal.II/base/config.h>

// The exception machinery (including the macros defined in
// exception_macros.h) references Kokkos functions. The places that
// use exceptions must know about these functions, and to avoid them
// all having to include Kokkos headers, we have to do it here:
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <Kokkos_Macros.hpp>
#if DEAL_II_KOKKOS_VERSION_GTE(4, 2, 0)
#  include <Kokkos_Abort.hpp>
#else
#  include <Kokkos_Core.hpp>
#endif
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS


#include <deal.II/base/exception_macros.h>

#include <complex>
#include <exception>
#include <ostream>
#include <string>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN


/**
 * This class is the base class for all exception classes. Do not use its
 * methods and variables directly since the interface and mechanism may be
 * subject to change. Rather create new exception classes using the
 * <tt>DeclException</tt> macro family.
 *
 * See the
 * @ref Exceptions
 * topic for more details on this class and what can be done with classes
 * derived from it.
 *
 * @ingroup Exceptions
 */
class ExceptionBase : public std::exception
{
public:
  /**
   * Default constructor.
   */
  ExceptionBase();

  /**
   * Copy constructor.
   */
  ExceptionBase(const ExceptionBase &exc);

  /**
   * Destructor.
   */
  virtual ~ExceptionBase() noexcept override = default;

  /**
   * Copy operator. This operator is deleted since exception objects
   * are not copyable.
   */
  ExceptionBase
  operator=(const ExceptionBase &) = delete;

  /**
   * Set the file name and line of where the exception appeared as well as the
   * violated condition and the name of the exception as a char pointer. This
   * function also populates the stacktrace.
   */
  void
  set_fields(const char *file,
             const int   line,
             const char *function,
             const char *cond,
             const char *exc_name);


  /**
   * Override the standard function that returns the description of the error.
   */
  virtual const char *
  what() const noexcept override;

  /**
   * Get exception name.
   */
  const char *
  get_exc_name() const;

  /**
   * Print out the general part of the error information.
   */
  void
  print_exc_data(std::ostream &out) const;

  /**
   * Print more specific information about the exception which occurred.
   * Overload this function in your own exception classes.
   */
  virtual void
  print_info(std::ostream &out) const;

  /**
   * Print a stacktrace, if one has been recorded previously, to the given
   * stream.
   */
  void
  print_stack_trace(std::ostream &out) const;

protected:
  /**
   * Name of the file this exception happens in.
   */
  const char *file;

  /**
   * Line number in this file.
   */
  unsigned int line;

  /**
   * Name of the function, pretty printed.
   */
  const char *function;

  /**
   * The violated condition, as a string.
   */
  const char *cond;

  /**
   * Name of the exception and call sequence.
   */
  const char *exc;

  /**
   * The number of stacktrace frames that are stored in the following variable.
   * Zero if the system does not support stack traces.
   */
  int n_stacktrace_frames;

#ifdef DEAL_II_HAVE_GLIBC_STACKTRACE
  /**
   * Array of pointers that contains the raw stack trace.
   */
  void *raw_stacktrace[25];
#endif

private:
  /**
   * Internal function that generates the c_string. Called by what().
   */
  void
  generate_message() const;

  /**
   * A pointer to the c_string that will be printed by what(). It is populated
   * by generate_message()
   */
  mutable std::string what_str;
};


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
   * @{
   */

  /**
   * Exception denoting a division by zero.
   */
  DeclExceptionMsg(ExcDivideByZero,
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
  DeclException1(
    ExcNumberNotFinite,
    std::complex<double>,
    << "In a significant number of places, deal.II checks that some intermediate "
    << "value is a finite number (as opposed to plus or minus infinity, or "
    << "NaN/Not a Number). In the current function, we encountered a number "
    << "that is not finite (its value is " << arg1 << " and therefore "
    << "violates the current assertion).\n\n"
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
    << "program visited before the place where this error is produced. "
    << "One way to check for finiteness is to use the 'AssertIsFinite' "
    << "macro.");

  /**
   * Trying to allocate a new object failed due to lack of free memory.
   */
  DeclException1(ExcOutOfMemory,
                 std::size_t,
                 "Your program tried to allocate some memory but this "
                 "allocation failed. Typically, this either means that "
                 "you simply do not have enough memory in your system, "
                 "or that you are (erroneously) trying to allocate "
                 "a chunk of memory that is simply beyond all reasonable "
                 "size, for example because the size of the object has "
                 "been computed incorrectly."
                 "\n\n"
                 "In the current case, the request was for "
                   << arg1 << " bytes.");

  /**
   * A memory handler reached a point where all allocated objects should have
   * been released. Since this exception is thrown, some were still allocated.
   */
  DeclException1(ExcMemoryLeak,
                 int,
                 << "Destroying memory handler while " << arg1
                 << " objects are still allocated.");

  /**
   * An error occurred reading or writing a file.
   */
  DeclExceptionMsg(ExcIO,
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
   * The constructor takes a single argument of type <tt>std::string</tt> naming
   * the file.
   */
  DeclException1(ExcFileNotOpen,
                 std::string,
                 << "Could not open file " << arg1
                 << "."
                    "\n\n"
                    "If this happens during an operation that tries to read "
                    "data: you may be "
                    "trying to read from a file that doesn't exist or that is "
                    "not readable given its file permissions."
                    "\n\n"
                    "If this happens during an operation that tries to write "
                    "data: you may be trying to write to a file to which file "
                    "or directory permissions do not allow you to write. A "
                    "typical example is where you specify an output file in "
                    "a directory that does not exist.");

  /**
   * Exception denoting a part of the library or application program that has
   * not yet been implemented. In many cases, this only indicates that there
   * wasn't much need for something yet, not that this is difficult to
   * implement. It is therefore quite worth the effort to take a look at the
   * corresponding place and see whether it can be implemented without too
   * much effort.
   */
  DeclExceptionMsg(ExcNotImplemented,
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
   * algorithm, these exceptions will indicate to them if they violate
   * assumptions that are used later in the algorithm. Furthermore, it
   * sometimes happens that an algorithm does not work in very rare corner
   * cases. These cases will then be trapped sooner or later by the exception,
   * so that the algorithm can then be fixed for these cases as well.
   */
  DeclExceptionMsg(ExcInternalError,
                   "This exception -- which is used in many places in the "
                   "library -- usually indicates that some condition which "
                   "the author of the code thought must be satisfied at a "
                   "certain point in an algorithm, is not fulfilled. An "
                   "example would be that the first part of an algorithm "
                   "sorts elements of an array in ascending order, and "
                   "a second part of the algorithm later encounters an "
                   "element that is not larger than the previous one."
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
  DeclExceptionMsg(
    ExcPureFunctionCalled,
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
   * This exception is used if some user function is not provided.
   */
  DeclException1(ExcFunctionNotProvided,
                 std::string,
                 << "Please provide an implementation for the function \""
                 << arg1 << "\"");

  /**
   * This exception is used if some user function returns nonzero exit codes.
   */
  DeclException2(
    ExcFunctionNonzeroReturn,
    std::string,
    int,
    << "The function \"" << arg1 << "\" returned the nonzero value " << arg2
    << ", but the calling site expected the return value to be zero. "
       "This error often happens when the function in question is a 'callback', "
       "that is a user-provided function called from somewhere within deal.II "
       "or within an external library such as PETSc, Trilinos, SUNDIALS, etc., "
       "that expect these callbacks to indicate errors via nonzero return "
       "codes.");

  /**
   * This exception is used if some object is found uninitialized.
   */
  DeclException0(ExcNotInitialized);

  /**
   * The object is in a state not suitable for this operation.
   */
  DeclException0(ExcInvalidState);

  /**
   * This exception is raised if a functionality is not possible in the given
   * dimension. Mostly used to throw function calls in 1d.
   *
   * The constructor takes a single <tt>int</tt>, denoting the dimension.
   */
  DeclException1(ExcImpossibleInDim,
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
  DeclException2(ExcImpossibleInDimSpacedim,
                 int,
                 int,
                 << "You are trying to execute functionality that is "
                 << "impossible in dimensions <" << arg1 << ',' << arg2
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
  DeclException2(ExcDimensionMismatch,
                 std::size_t,
                 std::size_t,
                 << "Two sizes or dimensions were supposed to be equal, "
                 << "but aren't. They are " << arg1 << " and " << arg2 << '.');

  /**
   * This exception is raised whenever deal.II cannot convert between integer
   * types.
   */
  DeclException2(
    ExcInvalidIntegerConversion,
    long long,
    long long,
    << "Two integers should be equal to each other after a type conversion but "
    << "aren't. A typical cause of this problem is that the integral types "
    << "used by deal.II and an external library are different (e.g., one uses "
    << "32-bit integers and the other uses 64-bit integers). The integers are "
    << arg1 << " and " << arg2 << '.');

  /**
   * The first dimension should be either equal to the second or the third,
   * but it is neither.
   */
  DeclException3(ExcDimensionMismatch2,
                 std::size_t,
                 std::size_t,
                 std::size_t,
                 << "The size or dimension of one object, " << arg1
                 << " was supposed to be "
                 << "equal to one of two values, but isn't. The two possible "
                 << "values are " << arg2 << " and " << arg3 << '.');

  /**
   * This exception indicates that an index is not within the expected range.
   * For example, it may be that you are trying to access an element of a
   * vector which does not exist.
   *
   * The constructor takes three <tt>std::size_t</tt> arguments, namely
   * <ol>
   * <li> the violating index
   * <li> the lower bound
   * <li> the upper bound plus one
   * </ol>
   */
  DeclException3(
    ExcIndexRange,
    std::size_t,
    std::size_t,
    std::size_t,
    << "Index " << arg1 << " is not in the half-open range [" << arg2 << ','
    << arg3 << ")."
    << (arg2 == arg3 ?
          " In the current case, this half-open range is in fact empty, "
          "suggesting that you are accessing an element of an empty "
          "collection such as a vector that has not been set to the "
          "correct size." :
          ""));

  /**
   * This exception indicates that an index is not within the expected range.
   * For example, it may be that you are trying to access an element of a
   * vector which does not exist.
   *
   * The constructor takes three arguments, namely
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
  DeclException3(
    ExcIndexRangeType,
    T,
    T,
    T,
    << "Index " << arg1 << " is not in the half-open range [" << arg2 << ','
    << arg3 << ")."
    << (arg2 == arg3 ?
          " In the current case, this half-open range is in fact empty, "
          "suggesting that you are accessing an element of an empty "
          "collection such as a vector that has not been set to the "
          "correct size." :
          ""));

  /**
   * A number is too small.
   */
  DeclException2(ExcLowerRange,
                 int,
                 int,
                 << "Number " << arg1 << " must be larger than or equal "
                 << arg2 << '.');

  /**
   * A generic exception definition for the ExcLowerRange above.
   */
  template <typename T>
  DeclException2(ExcLowerRangeType,
                 T,
                 T,
                 << "Number " << arg1 << " must be larger than or equal "
                 << arg2 << '.');

  /**
   * This exception indicates that the first argument should be an integer
   * multiple of the second, but is not.
   */
  DeclException2(ExcNotMultiple,
                 int,
                 int,
                 << "Division " << arg1 << " by " << arg2
                 << " has remainder different from zero.");

  /**
   * This exception is thrown if the iterator you access has corrupted data.
   * It might for instance be, that the container it refers does not have an
   * entry at the point the iterator refers.
   *
   * Typically, this will be an internal error of deal.II, because the
   * increment and decrement operators should never yield an invalid iterator.
   */
  DeclExceptionMsg(ExcInvalidIterator,
                   "You are trying to use an iterator, but the iterator is "
                   "in an invalid state. This may indicate that the iterator "
                   "object has not been initialized, or that it has been "
                   "moved beyond the end of the range of valid elements.");

  /**
   * This exception is thrown if the iterator you incremented or decremented
   * was already at its final state.
   */
  DeclExceptionMsg(ExcIteratorPastEnd,
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
  DeclException1(ExcMessage, std::string, << arg1);

  /**
   * Parallel vectors with ghost elements are read-only vectors.
   */
  DeclExceptionMsg(ExcGhostsPresent,
                   "You are trying an operation on a vector that is only "
                   "allowed if the vector has no ghost elements, but the "
                   "vector you are operating on does have ghost elements."
                   "\n\n"
                   "Specifically, there are two kinds of operations that "
                   "are typically not allowed on vectors with ghost elements. "
                   "First, vectors with ghost elements are read-only "
                   "and cannot appear in operations that write into these "
                   "vectors. Second, reduction operations (such as computing "
                   "the norm of a vector, or taking dot products between "
                   "vectors) are not allowed to ensure that each vector "
                   "element is counted only once (as opposed to once for "
                   "the owner of the element plus once for each process "
                   "on which the element is stored as a ghost copy)."
                   "\n\n"
                   "See the glossary entry on 'Ghosted vectors' for more "
                   "information.");

  /**
   * Exception indicating that one of the cells in the input to
   * Triangulation::create_triangulation() or a related function cannot be used.
   */
  DeclException1(ExcGridHasInvalidCell,
                 int,
                 << "Something went wrong when making cell " << arg1
                 << ". Read the docs and the source code "
                 << "for more information.");

  /**
   * Some of our numerical classes allow for setting all entries to zero using
   * the assignment operator <tt>=</tt>.
   *
   * In many cases, this assignment operator makes sense <b>only</b> for the
   * argument zero. In other cases, this exception is thrown.
   */
  DeclExceptionMsg(ExcScalarAssignmentOnlyForZeroValue,
                   "You are trying an operation of the form 'vector = C', "
                   "'matrix = C', or 'tensor = C' with a nonzero scalar value "
                   "'C'. However, such assignments are only allowed if the "
                   "C is zero, since the semantics for assigning any other "
                   "value are not clear. For example: one could interpret "
                   "assigning a matrix a value of 1 to mean the matrix has a "
                   "norm of 1, the matrix is the identity matrix, or the "
                   "matrix contains only 1s. Similar problems exist with "
                   "vectors and tensors. Hence, to avoid this ambiguity, such "
                   "assignments are not permitted.");

  /**
   * This function requires support for the LAPACK library.
   */
  DeclExceptionMsg(
    ExcNeedsLAPACK,
    "You are attempting to use functionality that is only available "
    "if deal.II was configured to use LAPACK, but when you configured "
    "the library, cmake did not find a valid LAPACK library."
    "\n\n"
    "You will have to ensure that your system has a usable LAPACK "
    "installation and re-install deal.II, making sure that cmake "
    "finds the LAPACK installation. You can check this by "
    "looking at the summary printed at the end of the cmake "
    "output.");

  /**
   * This function requires support for the HDF5 library.
   */
  DeclExceptionMsg(
    ExcNeedsHDF5,
    "You are attempting to use functionality that requires that deal.II is configured "
    "with HDF5 support. However, when you called 'cmake', HDF5 support "
    "was not detected."
    "\n\n"
    "You will have to ensure that your system has a usable HDF5 "
    "installation and re-install deal.II, making sure that cmake "
    "finds the HDF5 installation. You can check this by "
    "looking at the summary printed at the end of the cmake "
    "output.");

  /**
   * This function requires support for the MPI library.
   */
  DeclExceptionMsg(
    ExcNeedsMPI,
    "You are attempting to use functionality that is only available "
    "if deal.II was configured to use MPI."
    "\n\n"
    "You will have to ensure that your system has a usable MPI "
    "installation and re-install deal.II, making sure that cmake "
    "finds the MPI installation. You can check this by "
    "looking at the summary printed at the end of the cmake "
    "output.");

  /**
   * This function requires support for the FunctionParser library.
   */
  DeclExceptionMsg(
    ExcNeedsFunctionparser,
    "You are attempting to use functionality that is only available "
    "if deal.II was configured to use the function parser which "
    "relies on the muparser library, but cmake did not "
    "find a valid muparser library on your system and also did "
    "not choose the one that comes bundled with deal.II."
    "\n\n"
    "You will have to ensure that your system has a usable muparser "
    "installation and re-install deal.II, making sure that cmake "
    "finds the muparser installation. You can check this by "
    "looking at the summary printed at the end of the cmake "
    "output.");


  /**
   * This function requires support for the Assimp library.
   */
  DeclExceptionMsg(
    ExcNeedsAssimp,
    "You are attempting to use functionality that is only available "
    "if deal.II was configured to use Assimp, but cmake did not "
    "find a valid Assimp library."
    "\n\n"
    "You will have to ensure that your system has a usable Assimp "
    "installation and re-install deal.II, making sure that cmake "
    "finds the Assimp installation. You can check this by "
    "looking at the summary printed at the end of the cmake "
    "output.");

  /**
   * This function requires support for the Exodus II library.
   */
  DeclExceptionMsg(
    ExcNeedsExodusII,
    "You are attempting to use functionality that is only available if deal.II "
    "was configured to use Trilinos' SEACAS library (which provides ExodusII), "
    "but cmake did not find a valid SEACAS library."
    "\n\n"
    "You will have to ensure that your system has a usable ExodusII "
    "installation and re-install deal.II, making sure that cmake "
    "finds the ExodusII installation. You can check this by "
    "looking at the summary printed at the end of the cmake "
    "output.");

  /**
   * This function requires support for the CGAL library.
   */
  DeclExceptionMsg(
    ExcNeedsCGAL,
    "You are attempting to use functionality that is only available "
    "if deal.II was configured to use CGAL, but cmake did not "
    "find a valid CGAL library."
    "\n\n"
    "You will have to ensure that your system has a usable CGAL "
    "installation and re-install deal.II, making sure that cmake "
    "finds the CGAL installation. You can check this by "
    "looking at the summary printed at the end of the cmake "
    "output.");

#ifdef DEAL_II_WITH_MPI
  /**
   * Exception for MPI errors. This exception is only defined if
   * <code>deal.II</code> is compiled with MPI support. This exception should
   * be used with <code>AssertThrow</code> to check error codes of MPI
   * functions. For example:
   * @code
   * const int ierr = MPI_Isend(...);
   * AssertThrow(ierr == MPI_SUCCESS, ExcMPI(ierr));
   * @endcode
   * or, using the convenience macro <code>AssertThrowMPI</code>,
   * @code
   * const int ierr = MPI_Irecv(...);
   * AssertThrowMPI(ierr);
   * @endcode
   *
   * If the assertion fails then the error code will be used to print a helpful
   * message to the screen by utilizing the <code>MPI_Error_string</code>
   * function.
   *
   * @ingroup Exceptions
   */
  class ExcMPI : public dealii::ExceptionBase
  {
  public:
    ExcMPI(const int error_code);

    virtual void
    print_info(std::ostream &out) const override;

    const int error_code;
  };
#endif // DEAL_II_WITH_MPI



#ifdef DEAL_II_TRILINOS_WITH_SEACAS
  /**
   * Exception for ExodusII errors. This exception is only defined if
   * <code>deal.II</code> is compiled with SEACAS support, which is available
   * through Trilinos. This function should be used with the convenience macro
   * AssertThrowExodusII.
   *
   * @ingroup Exceptions
   */
  class ExcExodusII : public ExceptionBase
  {
  public:
    /**
     * Constructor.
     *
     * @param error_code The error code returned by an ExodusII function.
     */
    ExcExodusII(const int error_code);

    /**
     * Print a description of the error to the given stream.
     */
    virtual void
    print_info(std::ostream &out) const override;

    /**
     * Store the error code.
     */
    const int error_code;
  };
#endif // DEAL_II_TRILINOS_WITH_SEACAS

  /**
   * An exception to be thrown in user call-backs. See the glossary entry
   * on user call-back functions for more information.
   */
  DeclExceptionMsg(
    RecoverableUserCallbackError,
    "A user call-back function encountered a recoverable error, "
    "but the underlying library that called the call-back did not "
    "manage to recover from the error and aborted its operation."
    "\n\n"
    "See the glossary entry on user call-back functions for more "
    "information.");

  /** @} */

} /*namespace StandardExceptions*/



/**
 * In this namespace, functions in connection with the Assert and AssertThrow
 * mechanism are declared.
 *
 * @ingroup Exceptions
 */
namespace deal_II_exceptions
{
  namespace internals
  {
    /**
     * Setting this variable to false will disable deal.II's exception mechanism
     * to abort the problem. The Assert() macro will throw the exception instead
     * and the AssertNothrow() macro will just print the error message. This
     * variable should not be changed directly. Use disable_abort_on_exception()
     * instead.
     */
    extern bool allow_abort_on_exception;
  } // namespace internals

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
   *
   * @see Exceptions
   */
  void
  set_additional_assert_output(const char *const p);

  /**
   * Calling this function disables printing a stacktrace along with the other
   * output printed when an exception occurs. Most of the time, you will want
   * to see such a stacktrace; suppressing it, however, is useful if one wants
   * to compare the output of a program across different machines and systems,
   * since the stacktrace shows memory addresses and library names/paths that
   * depend on the exact setup of a machine.
   *
   * @see Exceptions
   */
  void
  suppress_stacktrace_in_exceptions();

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
   *
   * @see enable_abort_on_exception
   * @see Exceptions
   */
  void
  disable_abort_on_exception();

  /**
   * Calling this function switches on the use of <tt>std::abort()</tt> when
   * an exception is created using the Assert() macro, instead of throwing it.
   * This restores the standard behavior.
   *
   * @see disable_abort_on_exception
   * @see Exceptions
   */
  void
  enable_abort_on_exception();

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
     * Abort the program by printing the
     * error message provided by @p exc and calling <tt>std::abort()</tt>.
     */
    [[noreturn]] void
    abort(const ExceptionBase &exc) noexcept;

    /**
     * An enum describing how to treat an exception in issue_error_noreturn.
     */
    enum class ExceptionHandling
    {
      /**
       * Abort the program by calling <code>std::abort</code> unless
       * deal_II_exceptions::disable_abort_on_exception has been called: in
       * that case the program will throw an exception.
       */
      abort_or_throw_on_exception,
      /**
       * Throw the exception normally.
       */
      throw_on_exception
    };

    /**
     * This routine does the main work for the exception generation mechanism
     * used in the <tt>Assert</tt> and <tt>AssertThrow</tt> macros: as the
     * name implies, this function either ends by throwing an exception (if
     * @p handling is ExceptionHandling::throw_on_exception, or @p handling is try_abort_exception
     * and deal_II_exceptions::disable_abort_on_exception is false) or with a
     * call to <tt>abort</tt> (if @p handling is try_abort_exception and
     * deal_II_exceptions::disable_abort_on_exception is true).
     *
     * The actual exception object (the last argument) is typically an unnamed
     * object created in place; because we modify it, we can't take it by
     * const reference, and temporaries don't bind to non-const references.
     * So take it by value (=copy it) with a templated type to avoid slicing
     * -- the performance implications are pretty minimal anyway.
     *
     * @ref ExceptionBase
     */
    template <typename ExceptionType>
    [[noreturn]] void
    issue_error_noreturn(ExceptionHandling handling,
                         const char       *file,
                         int               line,
                         const char       *function,
                         const char       *cond,
                         const char       *exc_name,
                         ExceptionType     e)
    {
      static_assert(std::is_base_of_v<ExceptionBase, ExceptionType>,
                    "The provided exception must inherit from ExceptionBase.");
      // Fill the fields of the exception object
      e.set_fields(file, line, function, cond, exc_name);

      switch (handling)
        {
          case ExceptionHandling::abort_or_throw_on_exception:
            {
              if (dealii::deal_II_exceptions::internals::
                    allow_abort_on_exception)
                internals::abort(e);
              else
                {
                  // We are not allowed to abort, so just throw the error:
                  throw e;
                }
            }
          case ExceptionHandling::throw_on_exception:
            throw e;
          // this function should never return (and AssertNothrow can);
          // something must have gone wrong in the error handling code for us
          // to get this far, so throw an exception.
          default:
            throw ::dealii::StandardExceptions::ExcInternalError();
        }
    }

    /**
     * Internal function that does the work of issue_error_nothrow.
     */
    void
    do_issue_error_nothrow(const ExceptionBase &e) noexcept;

    /**
     * Exception generation mechanism in case we must not throw.
     *
     * @ref ExceptionBase
     *
     * @note This function is defined with a template for the same reasons as
     * issue_error_noreturn().
     */
    template <typename ExceptionType>
    void
    issue_error_nothrow(const char   *file,
                        int           line,
                        const char   *function,
                        const char   *cond,
                        const char   *exc_name,
                        ExceptionType e) noexcept
    {
      static_assert(std::is_base_of_v<ExceptionBase, ExceptionType>,
                    "The provided exception must inherit from ExceptionBase.");
      // Fill the fields of the exception object
      e.set_fields(file, line, function, cond, exc_name);
      // avoid moving a bunch of code into the header by dispatching to
      // another function:
      do_issue_error_nothrow(e);
    }
  } /*namespace internals*/

} /*namespace deal_II_exceptions*/



namespace deal_II_exceptions
{
  namespace internals
  {
    /**
     * A function that compares two values for equality, after converting to a
     * common type to avoid compiler warnings when comparing objects of
     * different types (e.g., unsigned and signed variables).
     */
    template <typename T, typename U>
    inline DEAL_II_HOST_DEVICE constexpr bool
    compare_for_equality(const T &t, const U &u)
    {
      using common_type = std::common_type_t<T, U>;
      return static_cast<common_type>(t) == static_cast<common_type>(u);
    }


    /**
     * A function that compares two values with `operator<`, after converting to
     * a common type to avoid compiler warnings when comparing objects of
     * different types (e.g., unsigned and signed variables).
     */
    template <typename T, typename U>
    inline DEAL_II_HOST_DEVICE constexpr bool
    compare_less_than(const T &t, const U &u)
    {
      using common_type = std::common_type_t<T, U>;
      return (static_cast<common_type>(t) < static_cast<common_type>(u));
    }
  } // namespace internals
} // namespace deal_II_exceptions



namespace internal
{
  // Workaround to allow for commas in template parameter lists
  // in preprocessor macros as found in
  // https://stackoverflow.com/questions/13842468/comma-in-c-c-macro
  template <typename F>
  struct argument_type;

  template <typename T, typename U>
  struct argument_type<T(U)>
  {
    using type = U;
  };

  template <typename F>
  using argument_type_t = typename argument_type<F>::type;
} // namespace internal


using namespace StandardExceptions;

DEAL_II_NAMESPACE_CLOSE

#endif
