// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


/**
 * @defgroup Exceptions Exceptions and assertions
 *
 * This module contains classes that are used in the exception mechanism of
 * deal.II.
 *
 * <h2>Brief overview</h2>
 * 
 * Exceptions are used in two different ways:
 * <ul>
 * 
 *   <li> Static assertions: These are checks that are only enabled in debug
 *   mode, not in optimized (or production) mode. They are meant to check that
 *   parameters to functions satisfy certain properties and similar
 *   assertions. For example, static assertions are used to make sure that two
 *   vectors that are added together have the same number of components --
 *   everything else would not make any sense anyway.
 *
 *   Such checks are performed by the Assert macro in several thousand places
 *   within the library. Also, several tutorial programs starting with step-5
 *   show how to do this.
 *
 *   If a static assertion is violated, the exception mechanism generates an
 *   exception of a type that indicates what exactly goes wrong, displays
 *   appropriate information, and then aborts the program -- if you try to add
 *   two vectors of different length, there is nothing that can be done within
 *   the program to cope with the situation, you have to go fix the program
 *   code instead. The exceptions of this module are used to indicate the
 *   reason for the failure.
 *
 *
 *   <li> Dynamic assertions: These are used to check dynamic features, such
 *   as whether an output file can be written to. These are things that can't
 *   be checked statically, i.e. they may change from program run to program
 *   run. It is therefore insufficient to only check these situations in debug
 *   mode.
 *
 *   Rather, one has to check them every time during execution of a
 *   program. Within deal.II, this is done using the AssertThrow macro
 *   introduced in step-9, step-13, and
 *   following tutorial programs. The macro checks a condition, and if
 *   violated throws an exception of one of the types declared in this
 *   module, using the C++ <code>throw</code> mechanism. Since these
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
 *  which by macro expansion does essentially the following:
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
 *  <h4>How it works internally</h4>
 *
 *  If the <tt>DEBUG</tt> preprocessor directive is set, the call <tt>Assert
 *  (cond, exc);</tt> is basically converted by the preprocessor into the
 *  following sequence:
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
 *  
 *  (Note that function names and exact calling sequences may change
 *  over time, but the general principle remains the same.) I.e., if
 *  the given condition is violated, then the file and line in which
 *  the exception occurred as well as the condition itself and the call
 *  sequence of the exception object is passed to the
 *  deal_II_exceptions::internals::issue_error_assert_1()
 *  function. Additionally an object of the form given by <tt>exc</tt>
 *  is created (this is normally an unnamed object like in
 *  <tt>ExcDomain (n, dim)</tt> of class <tt>ExcDomain</tt>) and
 *  transferred to this function.
 *
 *  <tt>__PRETTY__FUNCTION__</tt> is a macro defined by some compilers and
 *  gives the name of the function. If another compiler is used, we
 *  try to set this function to something reasonable, if the compiler
 *  provides us with that, and <tt>"(not available)"</tt> otherwise.
 *
 *  In <tt>issue_error_assert</tt>, the given data is transferred into
 *  the <tt>exc</tt> object by calling the set_fields() function;
 *  after that, the general error info is printed onto
 *  <tt>std::cerr</tt> using the PrintError() function of <tt>exc</tt>
 *  and finally the exception specific data is printed using the user
 *  defined function PrintError() (which is normally created using the
 *  <tt>DeclException (...)</tt> macro family. If it can be obtained
 *  from the operating system, the output may also contain a
 *  stacktrace to show where the error happened. Several of the
 *  @ref Tutorial programs show a typical output.
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
 *  <h3>Use of run-time exceptions</h3>
 *
 *  For this mode, the standard <tt>C++</tt> <tt>throw</tt> and <tt>catch</tt>
 *  concept exists. We
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
 *      std::cerr << "Exception occurred:" << std::endl
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
 *        e.set_fields (__FILE__, __LINE__, __PRETTY_FUNCTION__,
 *                      "condition as a string",
 *                      "name of condition as a string");
 *        throw e;
 *      };
 *  @endcode
 *
 *  For this purpose, the macro <tt>AssertThrow</tt> was invented. It does
 *  mainly the same job as does the <tt>Assert</tt> macro, but it does not
 *  kill the program, it rather throws an exception as shown above. The mode
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
 *  class name : public ExceptionBase {
 *    public:
 *      name (const type1 a1, const type2 a2) :
 *                     arg1 (a1), arg2(a2) {};
 *      virtual void print_info (std::ostream &out) const {
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
 * @author Wolfgang Bangerth, 1998-2006
 */
