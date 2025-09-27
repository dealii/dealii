// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * @defgroup Exceptions Exceptions and assertions
 *
 * This group contains classes that are used in the exception mechanism of
 * deal.II.
 *
 * <h2>Brief overview</h2>
 *
 * Exceptions are used in two different ways:
 * <ul>
 *
 *   <li> Static assertions: These are checks that are only enabled in debug
 *   mode, not in release (or optimized, production) mode. In deal.II, static
 *   assertions are typically used to check that parameters to functions satisfy
 *   certain properties, that internal data structures are consistent, and similar
 *   assertions. For example, static assertions are used to make sure that two
 *   vectors that are added together have the same number of components --
 *   everything else would not make any sense anyway.
 *
 *   Such checks are performed by the @p Assert macro in several thousand places
 *   within the library. Also, several tutorial programs starting with step-5
 *   show how to do this.
 *
 *   If a static assertion is violated, the exception mechanism generates an
 *   exception of a type that indicates what exactly goes wrong, displays
 *   appropriate information including the exact location where the problem
 *   was detected, and then aborts the program -- if you try to add
 *   two vectors of different length, there is nothing that can be done within
 *   the program to cope with the situation, you have to go fix the program
 *   code instead. There is generally not even a reason to @p throw an exception
 *   object using the usual C++ exception mechanism because there is nothing
 *   a function higher up could do in such cases to rectify the situation
 *   and deal with it in a useful way -- it's not that the program received
 *   bad data; the program is just buggy, and one can not intelligently
 *   work around that.
 *
 *   (It is sometimes useful to change the behavior of the @p Assert macro
 *   from aborting a program to throwing exceptions. On the other hand,
 *   exceptions are not allowed to propagate out of destructors of classes.
 *   For this purpose, there is a variant of the macro, called @p AssertNothrow,
 *   that can be used in destructors. These use cases are discussed further
 *   down below on this page.)
 *
 *
 *   <li> Dynamic assertions: These are used to check conditions that depend on
 *   external things that may be different from one program run to the next, such
 *   as whether an output file can be written to.
 *
 *   These are things that shouldn't
 *   be checked statically, because it is not guaranteed that a program for which
 *   the condition is satisfied in a debug mode run, will also have the condition
 *   satisfied in a subsequent release mode run -- in other words, it is not
 *   sufficient to only check these situations in debug mode.
 *
 *   Rather, one has to check them every time during execution of a
 *   program. Within deal.II, this is done using the @p AssertThrow macro
 *   introduced in step-9, step-13, and
 *   following tutorial programs. The macro checks a condition, and if
 *   violated throws an exception of one of the types declared in this
 *   group, using the C++ <code>throw</code> mechanism. Since these
 *   are run-time exceptions, this gives the program the chance to
 *   catch the exception and, for example, write the output to a
 *   writable file instead.
 * </ul>
 *
 *
 * <h2>Detailed description</h2>
 *
 *  The error handling mechanism in <tt>deal.II</tt> is generally used in two ways.
 *  The first uses error checking in debug mode only and is useful for programs
 *  which are not fully tested. When the program shows no errors anymore, one may
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
 *  in additional to the <tt>C++</tt> standard's <tt>std::exception</tt> class.
 *  Such a class is declared by the following lines of code:
 *  @code
 *     DeclException2 (ExcDomain, int, int,
 *                     << "Index= " << arg1 << "Upper Bound= " << arg2);
 *  @endcode
 *
 *  This declares an exception class named <tt>ExcDomain</tt>, which
 *  has two variables as additional information (named <tt>arg1</tt>
 *  and <tt>arg2</tt> by default) and which outputs the given sequence
 *  (which is appended to an <tt>std::ostream</tt> variable's name,
 *  thus the weird syntax). There are other <tt>DeclExceptionN</tt>
 *  macros for exception classes with more or no parameters. By
 *  convention, the name of all exception classes starts with
 *  <tt>Exc...</tt> and most of them are declared locally to the class
 *  it is to be used in (a few very frequently found ones are also
 *  declared in the StandardExceptions namespace and are available
 *  everywhere). Declaring exceptions globally is possible but
 *  pollutes the global namespace, is less readable and most of the time
 *  unnecessary.
 *
 *  Since exception classes are declared the same way for both modes
 *  of error checking, it is possible to use an exception declared
 *  through the <tt>DeclExceptionN(...)</tt> macro family for both
 *  static as well as dynamic checks.
 *
 *
 *  <h3>Use of the debug mode exceptions (static checks)</h3>
 *
 *  To use the exception mechanism for debug mode error checking, write lines
 *  like the following in your source code:
 *  @code
 *    Assert (n<dim, ExcDomain(n,dim));
 *  @endcode
 *  which by macro expansion does essentially the following (though the actual
 *  code is slightly more complicated):
 *  @code
 *    #ifdef DEBUG
 *        if (!(cond))
 *          {
 *            // issue error of class ExcDomain(n,dim)
 *          }
 *    #else
 *        // do nothing
 *    #endif
 *  @endcode
 *  That is, it issues an error only if the preprocessor variable
 *  <tt>DEBUG</tt> is set and if the given condition (in this case
 *  <tt>n < dim</tt> is violated).
 *
 *  If the exception was declared using the <tt>DeclException0 (...)</tt>
 *  macro, i.e., without any additional parameters, its name has
 *  nonetheless to be given with parentheses:
 *  <tt>Assert (i>m, ExcSomewhat());</tt>
 *
 *  <h4>How it works internally</h4>
 *
 *  If the <tt>DEBUG</tt> preprocessor directive is set, the call <tt>Assert
 *  (cond, exc);</tt> is basically converted by the preprocessor into the
 *  following sequence:
 *  @code
 *    if (!(cond))
 *      deal_II_exceptions::internals::issue_error_noreturn
 *             (__FILE__,
 *              __LINE__,
 *              __PRETTY_FUNCTION__,
 *              #cond,
 *              #exc,
 *              exc);
 *  @endcode
 *
 *  (Note that function names and exact calling sequences may change
 *  over time, but the general principle remains the same.) I.e., if
 *  the given condition is violated, then the file and line in which
 *  the exception occurred as well as the condition itself and the call
 *  sequence of the exception object is passed to the
 *  deal_II_exceptions::internals::issue_error_noreturn()
 *  function. Additionally an object of the form given by <tt>exc</tt>
 *  is created (this is normally an unnamed object like in
 *  <tt>ExcDomain (n, dim)</tt> of class <tt>ExcDomain</tt>) and
 *  transferred to this function.
 *
 *  <tt>__PRETTY_FUNCTION__</tt> is a macro defined by some compilers and
 *  gives the name of the function. If another compiler is used, we
 *  try to set this function to something reasonable, if the compiler
 *  provides us with that, and <tt>"(not available)"</tt> otherwise.
 *
 *  In <tt>issue_error_noreturn</tt>, the given data is transferred into the
 *  <tt>exc</tt> object by calling the set_fields() function; Afterwards the
 *  program is either aborted (and information about the exception is printed
 *  to deallog) or the exception is thrown. The <tt>Assert</tt> macro does the
 *  first path (print and abort); <tt>AssertThrow</tt> does the second
 *  (throw). This behavior is consistent with the descriptions of static and
 *  dynamic assertions earlier in this document. If it can be obtained from
 *  the operating system, the output may also contain a stacktrace to show
 *  where the error happened. Several of the @ref Tutorial programs show a
 *  typical output.
 *
 *  If the preprocessor variable <tt>DEBUG</tt> is not set then the
 *  <tt>Assert</tt> macro is expanded to <tt>{}</tt>.
 *
 *  Sometimes, there is no useful condition for an exception other
 *  than that the program flow should not have reached a certain point,
 *  e.g. a <tt>default</tt> section of a <tt>switch</tt> statement. In this case,
 *  raise the exception by the following construct:
 *  @code
 *    Assert (false, ExcInternalError());
 *  @endcode
 *  See the step-7 and several other of the tutorial programs for
 *  a use of this construct.
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
 *  <h3>Use of run-time exceptions (dynamic checks)</h3>
 *
 *  C++ has a mechanism to indicate that something exceptional has
 *  happened: exceptions that can be triggered by <tt>throw</tt> statements
 *  and captured by <tt>catch</tt> clauses, see for example
 *  https://en.wikipedia.org/wiki/C%2B%2B#Exception_handling and
 *  https://www.cplusplus.com/doc/tutorial/exceptions/ .
 *
 *  At some fundamental level, a typical C++ exception is an object that
 *  is placed in some special place, and then the function exits the current
 *  scope (e.g., the current function) through an exceptional return path.
 *  This is often enough to tell what problem triggered the exception,
 *  but more frequently it would be nice if one had more information: for
 *  example, in which line of the code the problem happened, or what
 *  non-existent entry of a sparse matrix the code wanted to write into.
 *
 *  Dynamic assertions in deal.II therefore extend this mechanism a bit.
 *  Typically, one would raise an exception by code such as
 *  @code
 *    if (!(cond))
 *      throw ExcSomething();
 *  @endcode
 *  and catch it using the statement
 *  @code
 *    try
 *      {
 *        do_something ();
 *      }
 *    catch (std::exception &e)
 *      {
 *        std::cerr << "Exception occurred:" << std::endl
 *                  << e.what ()
 *                  << std::endl;
 *        do_something_to_receiver ();
 *      }
 *  @endcode
 *  <tt>std::exception</tt> is a standard <tt>C++</tt> class providing basic functionality for
 *  exceptions, such as the virtual function <tt>what()</tt> that returns some
 *  information on the exception itself. This information is useful if an
 *  exception can't be handled properly, in which case as precise a description
 *  as possible should be printed.
 *
 *  The problem here is that to get significant and useful information out
 *  of <tt>what()</tt>, it is necessary to overload this function in our exception
 *  class and call the <tt>throw</tt> operator with additional arguments to the
 *  exception class. The first thing, overloading the <tt>what</tt> function is
 *  done using the <tt>DeclExceptionN</tt> macros, but putting the right information,
 *  which is the same as explained above for the <tt>Assert</tt> expansion, requires
 *  some work if one would want to write it down each time:
 *  @code
 *    if (!(cond))
 *      {
 *        ExcSomething e(additional information);
 *        e.set_fields (__FILE__, __LINE__, __PRETTY_FUNCTION__,
 *                      "condition as a string",
 *                      "name of condition as a string");
 *        throw e;
 *      };
 *  @endcode
 *
 *  For this purpose, the macro <tt>AssertThrow</tt> was invented. It does
 *  mainly the same job as does the <tt>Assert</tt> macro, but it does not
 *  abort the program; rather, it throws an exception as shown above. The mode
 *  of usage is
 *  @code
 *    AssertThrow (cond, ExcSomething(additional information));
 *  @endcode
 *
 *  The condition to be checked is incorporated into the macro in order to
 *  allow passing the violated condition as a string. The expansion of the
 *  <tt>AssertThrow</tt> macro is not affected by the <tt>DEBUG</tt>
 *  preprocessor variable.
 *
 *
 *  <h3>Description of the DeclExceptionN macro family</h3>
 *
 *  There is a whole family of <tt>DeclExceptionX</tt> macros
 *  where <tt>X</tt> is to be replaced by the number of additional
 *  parameters (0 to 5 presently).
 *  These macros are used to declare exception classes in the following
 *  way:
 *  @code
 *    DeclException2 (ExcDomain,
 *                    int,
 *                    int,
 *                    << " i=" << arg1 << ", m=" << arg2);
 *  @endcode
 *  The first argument denotes the name of the exception class to be created.
 *  The next arguments are the types of the parameters (in this
 *  case there two types, corresponding to the <tt>X</tt> in
 *  <tt>DeclExceptionX</tt>) and finally the output
 *  sequence with which you can print additional information.
 *
 *  The syntax of the output sequence is a bit weird but gets
 *  clearer once you see how this macro is defined (again schematically, actual
 *  function names and definitions may change over time and be different):
 *  @code
 *  class name : public ExceptionBase
 *  {
 *  public:
 *    name (const type1 a1, const type2 a2) : arg1 (a1), arg2(a2)
 *    {}
 *
 *    virtual void print_info (std::ostream &out) const
 *    {
 *      out << "    " outsequence << std::endl;
 *    }
 *  private:
 *    type1 arg1;
 *    type2 arg2;
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
 *  @code
 *    --------------------------------------------------------
 *    An error occurred in line <301> of file <exc-test.cc>.
 *    The violated condition was:
 *      i<m
 *    The name and call sequence of the exception was:
 *      MyExc2(i,m)
 *    Additional Information:
 *      i=5, m=3
 *    --------------------------------------------------------
 *  @endcode
 *
 *  Obviously for the <tt>DeclException0(name)</tt> macro, no types and
 *  also no output sequence is allowed.
 *
 *
 *  <h3>A corner case of @p Assert: The @p AssertNothrow macro</h3>
 *
 *  The default implementation of the @p Assert macro, as discussed above,
 *  prints detailed information about what exactly went wrong to the
 *  screen and then aborts the program. Aborting the program is useful
 *  because it allows easily finding the place where something went
 *  wrong -- including all of the information how we got to that
 *  place -- by running the program in a debugger.
 *
 *  On the other hand, there are cases where aborting a program may be
 *  undesirable and one needs to exit in a somewhat more graceful
 *  way -- even if there is really not very much one can do in these
 *  cases to still produce a meaningful result. An example is if a
 *  deal.II program is run as one module in a bigger framework of
 *  software. Think, for example, of a case where a deal.II program
 *  computed the flow field that corresponds to a set of input
 *  variables provided by some optimization routine: if the optimizer
 *  on the outside provided a negative density as input (a condition
 *  one might want to check via @p Assert), then this
 *  clearly makes no sense, and the flow solver cannot produce a
 *  meaningful answer; but it should tell that to the optimizer nicely,
 *  rather than just aborting the entire process (optimizer and flow
 *  solver).
 *
 *  For this purpose, one can call
 *  deal_II_exceptions::disable_abort_on_exception() that switches
 *  what @p Assert does from aborting the program to essentially the
 *  same as @p AssertThrow does, namely using the C++ @p throw mechanism
 *  to raise an exception. This exception can then be caught at a higher
 *  level -- e.g., in the optimization routine that sits atop the flow
 *  solver, and that can then decide what it wants to do with the
 *  situation.
 *
 *  This is all nice and good, but C++ does not allow throwing exceptions
 *  inside the destructors of classes, or in a function that is currently
 *  being called from a destructor higher up in the call stack. To this
 *  end, there is a separate macro, @p AssertNothrow, that can be used in
 *  destructors: It acts just like @p Assert usually does -- in particular,
 *  it only checks the condition in debug mode -- but it is immune to the
 *  effect of deal_II_exceptions::disable_abort_on_exception(): It will
 *  only ever abort the program, and never throw an exception.
 */
