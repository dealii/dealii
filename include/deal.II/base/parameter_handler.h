// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_parameter_handler_h
#define dealii_parameter_handler_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/subscriptor.h>

#include <boost/archive/basic_archive.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/property_tree/ptree_serialization.hpp>
#include <boost/serialization/split_member.hpp>

#include <map>
#include <memory>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations for interfaces and friendship
#ifndef DOXYGEN
class LogStream;
class MultipleParameterLoop;
#endif

/**
 * The ParameterHandler class provides a standard interface to an input file
 * which provides at run-time for program parameters such as time step sizes,
 * geometries, right hand sides etc. The input for the program is given in
 * files, streams or strings in memory using text like
 *   @code
 *     set Time step size = 0.3
 *     set Geometry       = [0,1]x[0,3]
 *   @endcode
 * Input may be sorted into subsection trees in order to give the input a
 * logical structure, and input files may include other files.
 *
 * The ParameterHandler class is discussed in step-29,
 * step-33, and step-34.
 *
 * <h3>Declaring entries</h3>
 *
 * In order to use the facilities of a ParameterHandler object, one first has
 * to make known the different entries the input file may or may not contain.
 * This is done in the following way:
 *
 *   @code
 *     ...
 *     ParameterHandler prm;
 *     prm.declare_entry ("Time step size",
 *                        "0.2",
 *                        Patterns::Double(),
 *                        "Some documentation");
 *     prm.declare_entry ("Geometry",
 *                        "[0,1]x[0,1]",
 *                        Patterns::Anything());
 *     ...
 *   @endcode
 * Each entry is declared using the function declare_entry(). The first
 * parameter is the name of the entry (in short: the entry). The second is the
 * default answer to be taken in case the entry is not specified in the input
 * file. The third parameter is a regular expression which the input (and the
 * default answer) has to match.  Several such regular expressions are defined
 * in Patterns. This parameter can be omitted, in which case it will default
 * to Patterns::Anything, i.e. a pattern that matches every input string. The
 * fourth parameter can be used to document the intent or expected format of
 * an entry; its value is printed as a comment when writing all entries of a
 * ParameterHandler object using the print_parameters() function to allow for
 * easier understanding of a parameter file. It can be omitted as well, in
 * which case no such documentation will be printed.
 *
 * Entries may be located in subsections which form a kind of input tree. For
 * example input parameters for linear solver routines should be classified in
 * a subsection named <tt>Linear solver</tt> or any other suitable name. This
 * is accomplished in the following way:
 * @code
 * ...
 * LinEq eq;
 * eq.declare_parameters (prm);
 * ...
 *
 * void LinEq::declare_parameters (ParameterHandler &prm)
 * {
 *   prm.enter_subsection("Linear solver");
 *   {
 *     prm.declare_entry ("Solver",
 *                        "CG",
 *                        Patterns::Selection("CG|GMRES|GaussElim"),
 *                        "Name of a linear solver for the inner iteration");
 *     prm.declare_entry ("Maximum number of iterations", "20",
 *                        ParameterHandler::RegularExpressions::Integer());
 *     ...
 *   }
 *   prm.leave_subsection ();
 * }
 * @endcode
 *
 * Subsections may be nested. For example a nonlinear solver may have a linear
 * solver as member object. Then the function call tree would be something
 * like (if the class <tt>NonLinEq</tt> has a member variables <tt>eq</tt> of
 * type <tt>LinEq</tt>):
 * @code
 * void NonLinEq::declare_parameters (ParameterHandler &prm)
 * {
 *   prm.enter_subsection ("Nonlinear solver");
 *   {
 *     prm.declare_entry ("Nonlinear method",
 *                        "Newton-Raphson",
 *                        ParameterHandler::RegularExpressions::Anything());
 *     eq.declare_parameters (prm);
 *   }
 *   prm.leave_subsection ();
 * }
 * @endcode
 *
 * For class member functions which declare the different entries we propose
 * to use the common name <tt>declare_parameters</tt>. In normal cases this
 * method can be <tt>static</tt> since the entries will not depend on any
 * previous knowledge. Classes for which entries should logically be grouped
 * into subsections should declare these subsections themselves. If a class
 * has two or more member variables of the same type both of which should have
 * their own parameters, this parent class' method <tt>declare_parameters</tt>
 * is responsible to group them into different subsections:
 * @code
 * void NonLinEq::declare_parameters (ParameterHandler &prm)
 * {
 *   prm.enter_subsection ("Nonlinear solver");
 *   {
 *     prm.enter_subsection ("Linear solver 1");
 *     {
 *       eq1.declare_parameters (prm);
 *     }
 *     prm.leave_subsection ();
 *
 *     prm.enter_subsection ("Linear solver 2");
 *     {
 *       eq2.declare_parameters (prm);
 *     }
 *     prm.leave_subsection ();
 *   }
 *   prm.leave_subsection ();
 * }
 * @endcode
 *
 *
 * <h3>Input files and special characters</h3>
 *
 * For the first example above the input file would look like the following:
 *   @code
 *     ...
 *     subsection Nonlinear solver
 *       set Nonlinear method = Gradient
 *       # this is a comment
 *       subsection Linear solver
 *         set Solver                       = CG
 *         set Maximum number of iterations = 30
 *       end
 *     end
 *     ...                       # other stuff
 *   @endcode
 * The words <tt>subsection</tt>, <tt>set</tt> and <tt>end</tt> may be either
 * written in lowercase or uppercase letters. Leading and trailing whitespace
 * is removed, multiple whitespace is condensed into only one. Since the
 * latter applies also to the name of an entry, an entry name will not be
 * recognized if in the declaration multiple whitespace is used.
 *
 * In entry names and values the following characters are not allowed:
 * <tt>\#</tt>, <tt>{</tt>, <tt>}</tt>, <tt>|</tt>. Their use is reserved for
 * the MultipleParameterLoop class.
 *
 * Comments starting with \# are skipped.
 *
 * Continuation lines are allowed by means of the character <tt>\\</tt>, which
 * must be the last character (aside from whitespace, which is ignored) of the
 * line. When a line is a continuation (i.e., the previous line ended in a
 * <tt>\\</tt>), then, unlike the default behavior of the <tt>C</tt>
 * preprocessor, all whitespace at the beginning of the line is ignored.
 *
 * We propose to use the following scheme to name entries: start the first
 * word with a capital letter and use lowercase letters further on. The same
 * applies to the possible entry values to the right of the <tt>=</tt> sign.
 *
 *
 * <h3>Including other input files</h3>
 *
 * An input file can include other include files using the syntax
 *   @code
 *     ...
 *     include some_other_file.prm
 *     ...
 *   @endcode
 * The file so referenced is searched for relative to the current directory
 * (not relative to the directory in which the including parameter file is
 * located, since this is not known to all three versions of the parse_input()
 * function).
 *
 *
 * <h3>Reading data from input sources</h3>
 *
 * In order to read input there are three possibilities: reading from an
 * <tt>std::istream</tt> object, reading from a file of which the name is
 * given and reading from a string in memory in which the lines are separated
 * by <tt>@\n</tt> characters. These possibilities are used as follows:
 *   @code
 *     ParameterHandler prm;
 *     ...
 *     // declaration of entries
 *     ...
 *     prm.parse_input (std::cin); // read input from standard in,
 *     // or
 *     prm.parse_input ("simulation.prm");
 *     // or
 *     char *in = "set Time step size = 0.3 \n ...";
 *     prm.parse_input_from_string (in);
 *     ...
 *   @endcode
 * You can use several sources of input successively. Entries which are
 * changed more than once will be overwritten every time they are used.
 *
 * You should not try to declare entries using declare_entry() and
 * enter_subsection() with as yet unknown subsection names after using
 * parse_input(). The results in this case are unspecified.
 *
 * If an error occurs upon reading the input, error messages are written to
 * <tt>std::cerr</tt> and the reader function returns with a return value of
 * <code>false</code>. This is opposed to almost all other functions in
 * deal.II, which would normally throw an exception if an error occurs; this
 * difference in behavior is a relic of the fact that this class predates
 * deal.II and had previously been written for a different project.
 *
 *
 * <h3>Using the %ParameterHandler Graphical User Interface</h3>
 *
 * An alternative to using the hand-written input files shown above is to use
 * the graphical user interface (GUI) that accompanies this class.
 *
 * See <a href="https://github.com/dealii/parameter_gui">the parameter_gui
 * github repository</a> for further details.
 *
 * <h3>Getting entry values out of a %ParameterHandler object</h3>
 *
 * Each class gets its data out of a ParameterHandler object by calling the
 * get()  member functions like this:
 *   @code
 *     void NonLinEq::get_parameters (ParameterHandler &prm)
 *     {
 *       prm.enter_subsection ("Nonlinear solver");
 *       std::string method = prm.get ("Nonlinear method");
 *       eq.get_parameters (prm);
 *       prm.leave_subsection ();
 *     }
 *   @endcode
 * get() returns the value of the given entry. If the entry was not specified
 * in the input source(s), the default value is returned. You have to enter
 * and leave subsections exactly as you did when declaring subsections. You may
 * choose the order in which to traverse the subsection tree.
 *
 * It is possible to avoid calls to enter_subsection() and leave_subsection()
 * by supplying get() with a vector of strings representing the path from
 * which to get a value. For example, the following two versions of
 * get_parameters() will produce the same result:
 *   @code
 *     void NonLinEq::get_parameters (ParameterHandler &prm)
 *     {
 *       prm.enter_subsection ("Equation 1 Settings");
 *       prm.enter_subsection ("Linear solver");
 *       solver_ = prm.get ("Solver");
 *       prm.leave_subsection ();
 *       prm.leave_subsection ();
 *     }
 *   @endcode
 *
 *   @code
 *     void NonLinEq::get_parameters (const ParameterHandler &prm)
 *     {
 *       std::vector<std::string> path =
 *         {"Equation 1 Settings", "Linear solver"};
 *       solver_ = prm.get (path, "Solver");
 *     }
 *   @endcode
 *
 * The latter method allows the ParameterHandler reference to be @p const.
 *
 * It is guaranteed that only entries matching the given regular expression
 * are returned, i.e. an input entry value which does not match the regular
 * expression is not stored.
 *
 * You can use get() to retrieve the parameter in text form, get_integer() to
 * get an integer or get_double() to get a double. You can also use
 * get_bool(). It will cause an internal error if the string could not be
 * converted to an integer, double or a bool. This should, though, not happen
 * if you correctly specified the regular expression for this entry; you
 * should not try to get out an integer or a double from an entry for which no
 * according regular expression was set. The internal error is raised through
 * the Assert() macro family which only works in debug mode.
 *
 * If you want to print out all user selectable features, use the
 * print_parameters() function. It is generally a good idea to print all
 * parameters at the beginning of a log file, since this way input and output
 * are together in one file which makes matching at a later time easier.
 * Additionally, the function also print those entries which have not been
 * modified in the input file and are thus set to default values; since
 * default values may change in the process of program development, you cannot
 * know the values of parameters not specified in the input file.
 *
 *
 *
 * <h3>Adding Actions to Parameters</h3>
 *
 * It is often convenient to have something happen as soon as a parameter
 * value is read. This could be a check that it is valid -- say, that a
 * file that is listed in the parameter file exists -- or to initiate
 * something else in response, such as setting a variable outside the
 * ParameterHandler (as in the example shown below). In almost all cases,
 * this "action" could also be initiated once all parameters are read
 * via parse_input(), but it is sometimes <i>convenient</i> to do it
 * right away.
 *
 * This is facilitated by the add_action() function that can be called
 * after declaring a parameter via declare_entry(). "Actions" are in essence
 * pointers to functions that will be called for parameters that have
 * associated actions. These functions take the value of a parameter as
 * argument, and can then do whatever they want with it -- e.g., save it
 * somewhere outside the ParameterHandler object. (Exactly when the
 * action is called is described in the documentation of the
 * add_action() function.) Of course, in C++ one doesn't usually pass
 * around the address of a function, but an action can be a function-like
 * object (taking a string as argument) that results from calling
 * such as a
 * <a href="http://en.cppreference.com/w/cpp/language/lambda">lambda
 * function</a> that has the form
 * @code
 *   [] (const std::string &value) { ... do something with the value ... }
 * @endcode
 * and that is attached to a specific parameter.
 *
 * A typical example of such an action would be as follows: let's assume
 * that you have a program that declares a parameter for the number
 * of iterations it is going to run, say
 * @code
 *   class MyAlgorithm
 *   {
 *      public:
 *        void run ();
 *      private:
 *        unsigned int n_iterations;
 *   };
 * @endcode
 * then one could obtain this parameter from a parameter file using a
 * code snippet in @p run() as follows:
 * @code
 *   void MyAlgorithm::run ()
 *   {
 *     ParameterHandler prm;
 *     prm.declare_entry ("Number of iterations",  // name of parameter
 *                        "10",                    // default value
 *                        Patterns::Integer(1,100),// allowed values: 1...100
 *                        "The number of ...");    // some documentation
 *
 *     // next read the parameter from an input file...
 *     prm.parse_input ("my_algorithm.prm");
 *
 *     // ...and finally get the value for use in the program:
 *     n_iterations = prm.get_integer ("Number of iterations");
 *
 *     ... actual code doing something useful follows here...
 * @endcode
 *
 * This two-step process -- first declaring the parameter, and later reading
 * it -- is a bit cumbersome because one has to first declare <i>all</i>
 * parameters and at a later time retrieve them from the ParameterHandler
 * object. In large programs, these two things also often happen in
 * different functions.
 *
 * To avoid this, it would be nice if we could put both the declaration
 * and the retrieval into the same place. This can be done via actions,
 * and the function would then look like this:
 * @code
 *   void MyAlgorithm::run ()
 *   {
 *     ParameterHandler prm;
 *     prm.declare_entry ("Number of iterations",  // name of parameter
 *                        "10",                    // default value
 *                        Patterns::Integer(1,100),// allowed values: 1...100
 *                        "The number of ...");    // some documentation
 *     prm.add_action ("Number of iterations",
 *                     [&](const std::string &value)
 *                     {
 *                       this->n_iterations = Utilities::string_to_int(value);
 *                     });
 *
 *     // next read the parameter from an input file...
 *     prm.parse_input ("my_algorithm.prm");
 *
 *     ... actual code doing something useful follows here...
 * @endcode
 * Here, the action consists of a lambda function that takes the value
 * for this parameter as a string, and then converts it to an integer
 * to store in the variable where it belongs. This action is
 * executed inside the call to <code>prm.parse_input()</code>, and so
 * there is now no longer a need to extract the parameter's value
 * at a later time. Furthermore, the code that sets the member variable
 * is located right next to the place where the parameter is actually
 * declared, so we no longer need to have two separate parts of the code
 * base that deal with input parameters.
 *
 * Of course, it is possible to execute far more involved actions than
 * just setting a member variable as shown above, even though that is
 * a typical case.
 *
 *
 * <h3>Style guide for data retrieval</h3>
 *
 * We propose that every class which gets data out of a ParameterHandler
 * object provides a function named <tt>get_parameters</tt>. This should be
 * declared <tt>virtual</tt>. <tt>get_parameters</tt> functions in derived
 * classes should call the <tt>BaseClass::get_parameters</tt> function.
 *
 *
 * <h3>Experience with large parameter lists</h3>
 *
 * Experience has shown that in programs defining larger numbers of parameters
 * (more than, say, fifty) it is advantageous to define an additional class
 * holding these parameters. This class is more like a C-style structure,
 * having a large number of variables, usually public. It then has at least
 * two functions, which declare and parse the parameters. In the main program,
 * the main class has an object of this parameter class and delegates
 * declaration and parsing of parameters to this object.
 *
 * The advantage of this approach is that you can keep out the technical
 * details (declaration and parsing) out of the main class and additionally
 * don't clutter up your main class with dozens or more variables denoting the
 * parameters.
 *
 *
 *
 * <h3>Worked Example</h3>
 *
 * This is the code:
 *   @code
 *   #include <deal.II/base/parameter_handler.h>
 *
 *   #include <iostream>
 *   #include <string>
 *
 *   using namespace dealii;
 *   class LinearEquation
 *   {
 *   public:
 *     static void declare_parameters (ParameterHandler &prm);
 *     void get_parameters (ParameterHandler &prm);
 *   private:
 *     std::string method;
 *     int         max_iterations;
 *   };
 *
 *
 *
 *   class Problem
 *   {
 *   private:
 *     LinearEquation eq1, eq2;
 *     std::string matrix1, matrix2;
 *     std::string outfile;
 *   public:
 *     static void declare_parameters (ParameterHandler &prm);
 *     void get_parameters (ParameterHandler &prm);
 *
 *     void do_something ();
 *   };
 *
 *
 *
 *   void LinearEquation::declare_parameters (ParameterHandler &prm)
 *   {
 *     // declare parameters for the linear solver in a subsection
 *     prm.enter_subsection ("Linear solver");
 *     {
 *       prm.declare_entry ("Solver",
 *                          "CG",
 *                          Patterns::Selection("CG|BiCGStab|GMRES"),
 *                          "Name of a linear solver for the inner iteration");
 *       prm.declare_entry ("Maximum number of iterations",
 *                          "20",
 *                          Patterns::Integer());
 *     }
 *     prm.leave_subsection ();
 *   }
 *
 *
 *
 *   void LinearEquation::get_parameters (ParameterHandler &prm)
 *   {
 *     prm.enter_subsection ("Linear solver");
 *     {
 *       method         = prm.get ("Solver");
 *       max_iterations = prm.get_integer ("Maximum number of iterations");
 *     }
 *     prm.leave_subsection ();
 *     std::cout << "  LinearEquation: method=" << method
 *               << ", max_iterations=" << max_iterations
 *               << std::endl;
 *   }
 *
 *
 *
 *   void Problem::declare_parameters (ParameterHandler &prm)
 *   {
 *     // first some global parameter entries
 *     prm.declare_entry (
 *       "Output file",
 *       "out",
 *       Patterns::Anything(),
 *       "Name of the output file, either relative or absolute");
 *     prm.declare_entry ("Equation 1", "Laplace",
 *                        Patterns::Anything(),
 *                        "String identifying the equation we want to solve");
 *     prm.declare_entry ("Equation 2",
 *                        "Elasticity",
 *                        Patterns::Anything());
 *
 *     // declare parameters for the first equation
 *     prm.enter_subsection ("Equation 1 Settings");
 *     {
 *       prm.declare_entry ("Matrix type",
 *                          "Sparse",
 *                          Patterns::Selection("Full|Sparse|Diagonal"),
 *                          "Type of the matrix to be used, either full, "
 *                          "sparse, or diagonal");
 *       LinearEquation::declare_parameters (prm);  // for eq1
 *     }
 *     prm.leave_subsection ();
 *
 *     // declare parameters for the second equation
 *     prm.enter_subsection ("Equation 2 Settings");
 *     {
 *       prm.declare_entry ("Matrix type",
 *                          "Sparse",
 *                          Patterns::Selection("Full|Sparse|Diagonal"));
 *       LinearEquation::declare_parameters (prm);  // for eq2
 *     }
 *     prm.leave_subsection ();
 *   }
 *
 *
 *
 *   void Problem::get_parameters (ParameterHandler &prm)
 *   {
 *     // entries of the problem class
 *     outfile = prm.get ("Output file");
 *     std::string equation1 = prm.get ("Equation 1"),
 *                 equation2 = prm.get ("Equation 2");
 *
 *     // get parameters for the first equation
 *     prm.enter_subsection ("Equation 1 Settings");
 *     {
 *       matrix1 = prm.get ("Matrix type");
 *       eq1.get_parameters (prm); // for eq1
 *     }
 *     prm.leave_subsection ();
 *
 *     // get parameters for the second equation
 *     prm.enter_subsection ("Equation 2 Settings");
 *     {
 *       matrix2 = prm.get ("Matrix type");
 *       eq2.get_parameters (prm); // for eq2
 *     }
 *     prm.leave_subsection ();
 *     std::cout
 *       << "  Problem: outfile=" << outfile << '\n'
 *       << "           eq1="     << equation1 << ", eq2=" << equation2 << '\n'
 *       << "           matrix1=" << matrix1 << ", matrix2=" << matrix2
 *       << std::endl;
 *   }
 *
 *
 *
 *   void Problem::do_something ()
 *   {
 *     // While this example does nothing here, at this point in the program
 *     // all of the parameters are known so we can start doing computations.
 *   }
 *
 *
 *
 *   int main ()
 *   {
 *     ParameterHandler prm;
 *     Problem p;
 *     p.declare_parameters (prm);
 *     // read input from "prmtest.prm"; giving argv[1] would also be a
 *     // good idea
 *     prm.parse_input ("prmtest.prm");
 *     // print parameters to std::cout as ASCII text
 *     std::cout << "\n\n";
 *     prm.print_parameters (std::cout, ParameterHandler::Text);
 *     // get parameters into the program
 *     std::cout << "\n\n" << "Getting parameters:" << std::endl;
 *     p.get_parameters (prm);
 *     // now run the program with these input parameters
 *     p.do_something ();
 *   }
 *   @endcode
 *
 *
 * This is the input file (named "prmtest.prm"):
 *   @code
 *   # first declare the types of equations
 *   set Equation 1 = Poisson
 *   set Equation 2 = Stokes
 *
 *   subsection Equation 1 Settings
 *     set Matrix type = Sparse
 *     subsection Linear solver # parameters for linear solver 1
 *       set Solver                       = Gauss-Seidel
 *       set Maximum number of iterations = 40
 *     end
 *   end
 *
 *   subsection Equation 2 Settings
 *     set Matrix type = Full
 *     subsection Linear solver
 *       set Solver                       = CG
 *       set Maximum number of iterations = 100
 *     end
 *   end
 *   @endcode
 *
 * And here is the output of the program:
 *   @code
 *   Line <8> of file <prmtest.prm>:
 *       The entry value
 *           Gauss-Seidel
 *       for the entry named
 *           Solver
 *       does not match the given pattern
 *           [Selection CG|BiCGStab|GMRES ]
 *
 *
 *   # Listing of Parameters
 *   # ---------------------
 *   # String identifying the equation we want to solve
 *   set Equation 1  = Poisson # default: Laplace
 *   set Equation 2  = Stokes  # default: Elasticity
 *
 *   # Name of the output file, either relative to the present path or absolute
 *   set Output file = out
 *
 *
 *   subsection Equation 1 Settings
 *     # Type of the matrix to be used, either full, sparse, or diagonal
 *     set Matrix type = Sparse
 *
 *
 *     subsection Linear solver
 *       set Maximum number of iterations = 40 # default: 20
 *       # Name of a linear solver for the inner iteration
 *       set Solver                       = CG
 *     end
 *
 *   end
 *
 *
 *   subsection Equation 2 Settings
 *     set Matrix type = Full # default: Sparse
 *
 *
 *     subsection Linear solver
 *       set Maximum number of iterations = 100 # default: 20
 *       # Name of a linear solver for the inner iteration
 *       set Solver                       = CG
 *     end
 *
 *   end
 *
 *
 *
 *
 *   Getting parameters:
 *     LinearEquation: method=CG, max_iterations=40
 *     LinearEquation: method=CG, max_iterations=100
 *     Problem: outfile=out
 *              eq1=Poisson, eq2=Stokes
 *              matrix1=Sparse, matrix2=Full
 *   @endcode
 *
 *
 *
 * <h3>Representation of Parameters</h3>
 *
 * Here is some more internal information about the representation of
 * parameters:
 *
 * Logically, parameters and the nested sections they are arranged in can be
 * thought of as a hierarchical directory structure, or a tree. Take, for
 * example, the following code declaring a set of parameters and sections they
 * live in:
 *   @code
 *     ParameterHandler prm;
 *
 *     prm.declare_entry ("Maximal number of iterations",
 *                        "10",
 *                        Patterns::Integer (1, 1000),
 *                        "A parameter that describes the maximal number of "
 *                        "iterations the CG method is to take before giving "
 *                        "up on a matrix.");
 *     prm.enter_subsection ("Preconditioner");
 *     {
 *       prm.declare_entry(
 *         "Kind",
 *         "SSOR",
 *         Patterns::Selection ("SSOR|Jacobi"),
 *         "A string that describes the kind of preconditioner to use.");
 *       prm.declare_entry(
 *         "Relaxation factor",
 *         "1.0",
 *         Patterns::Double (0, 1),
 *         "The numerical value (between zero and one) for the "
 *         "relaxation factor to use in the preconditioner.");
 *     }
 *     prm.leave_subsection ();
 *   @endcode
 *
 * We can think of the parameters so arranged as a file system in which every
 * parameter is a directory. The name of this directory is the name of the
 * parameter, and in this directory lie files that describe the parameter.
 * These files are at the time of writing this documentation (other fields,
 * such as those indicating "actions" may also exist in each directory):
 *
 * - <code>value</code>: The content of this file is the current value of this
 * parameter; initially, the content of the file equals the default value of
 * the parameter.
 * - <code>default_value</code>: The content of this file is the default value
 * of the parameter.
 * - <code>pattern</code>: A textual representation of the pattern that
 * describes the parameter's possible values.
 * - <code>pattern_index</code>: A number that indexes the Patterns::PatternBase
 * object that is used to describe the parameter.
 * - <code>documentation</code>: The content of this file is the documentation
 * given for a parameter as the last argument of the
 * ParameterHandler::declare_entry call. With the exception of the
 * <code>value</code> file, the contents of files are never changed after
 * declaration of a parameter.
 *
 * Alternatively, a directory in this file system may not have a file called
 * <code>value</code> in it. In that case, the directory represents a
 * subsection as declared above, and the directory's name will correspond to
 * the name of the subsection. It will then have no files in it at all, but it
 * may have further directories in it: some of these directories will be
 * parameters (indicates by the presence of files) or further nested
 * subsections.
 *
 * Given this explanation, the code above will lead to a hierarchical
 * representation of data that looks like this (the content of files is
 * indicated at the right in a different font):
 *
 * @image html parameter_handler.png
 *
 * Once parameters have been read in, the contents of the <code>value</code>
 * "files" may be different while the other files remain untouched.
 *
 * Using the ParameterHandler::print_parameters() function with
 * ParameterHandler::XML as second argument, we can get a complete
 * representation of this data structure in XML. It will look like this:
 *   @code
 *   <?xml version="1.0" encoding="utf-8"?>
 *   <ParameterHandler>
 *     <Maximal_20number_20of_20iterations>
 *       <value>10</value>
 *       <default_value>10</default_value>
 *       <documentation>
 *         A parameter that describes the maximal number of iterations the CG
 *         method is to take before giving up on a matrix.
 *       </documentation>
 *       <pattern>0</pattern>
 *       <pattern_description>
 *         [Integer range 1...1000 (inclusive)]
 *       </pattern_description>
 *     </Maximal_20number_20of_20iterations>
 *     <Preconditioner>
 *       <Kind><value>SSOR</value>
 *         <default_value>SSOR</default_value>
 *         <documentation>
 *           A string that describes the kind of preconditioner to use.
 *         </documentation>
 *         <pattern>1</pattern>
 *         <pattern_description>SSOR|Jacobi</pattern_description>
 *       </Kind>
 *       <Relaxation_20factor>
 *         <value>1.0</value>
 *         <default_value>1.0</default_value>
 *         <documentation>
 *           The numerical value (between zero and one) for the relaxation
 *           factor to use in the preconditioner.
 *         </documentation>
 *         <pattern>2</pattern>
 *         <pattern_description>
 *           [Floating point range 0...1 (inclusive)]
 *         </pattern_description>
 *       </Relaxation_20factor>
 *     </Preconditioner>
 *   <ParameterHandler>
 *   @endcode
 * This representation closely resembles the directory/file structure
 * discussed above. The only difference is that directory and file names are
 * mangled: since they should only contain letters and numbers, every
 * character in their names that is not a letter or number is replaced by an
 * underscore followed by its two-digit hexadecimal representation. In
 * addition, the special name "value" is mangled when used as the name of a
 * parameter, given that this name is also used to name special files in the
 * hierarchy structure. Finally, the entire tree is wrapped into a tag
 * <code>%ParameterHandler</code> to satisfy the XML requirement that there be
 * only a single top-level construct in each file.
 *
 * The tree structure (and its XML representation) is what the graphical user
 * interface (see above) uses to represent parameters like a directory/file
 * collection.
 *
 *
 * @ingroup input
 * @author Wolfgang Bangerth, October 1997, revised February 1998, 2010, 2011, 2017
 * @author Alberto Sartori, 2015
 * @author David Wells, 2016
 * @author Denis Davydov, 2018
 * @author Pasquale Claudio Africa, 2020
 */
class ParameterHandler : public Subscriptor
{
public:
  /**
   * List of possible output formats used for functions like
   * ParameterHandler::print_parameters(). The options can be categorized into
   * two groups:
   * - format options: PRM, LaTeX, Description, XML, JSON
   * - stylistic options: Short, KeepDeclarationOrder
   *
   * Only one format option may be specified at the time. Any function that
   * accepts an OutputStyle as an option will throw if you specify more than
   * one.
   *
   * A number of shortcuts of commonly used option combinations are provided.
   * E.g., ShortPRM prints the parameters in the PRM format, while skipping the
   * documentation.
   */
  enum OutputStyle
  {
    /**
     * Default stylistic style: print documentation and sort all parameters
     * alphabetically.
     */
    DefaultStyle = 0x0000,

    /**
     * Write input for ParameterHandler without comments or changed default
     * values.
     */
    Short = 0x0001,

    /**
     * Keep the order of the parameters as they have been declared.
     */
    KeepDeclarationOrder = 0x0002,

    /**
     * Write human readable output suitable to be read by
     * ParameterHandler::parse_input() again.
     */
    PRM = 0x0010,

    /**
     * Write human readable output suitable to be read by
     * ParameterHandler::parse_input() again.
     *
     * @deprecated Use `PRM` instead of `Text`.
     */
    Text = PRM,

    /**
     * Write parameters as a LaTeX table.
     */
    LaTeX = 0x0020,

    /**
     * Write out declared parameters with description and possible values.
     *
     * @note This format is not suitable to be read back again.
     */
    Description = 0x0040,

    /**
     * Write out everything as an <a
     * href="http://en.wikipedia.org/wiki/XML">XML</a> file suitable to be read
     * by ParameterHandler::parse_input_from_xml() again.
     *
     * See the general documentation of this class for an example of output.
     */
    XML = 0x0080,

    /**
     * Write out everything as a <a
     * href="http://en.wikipedia.org/wiki/JSON">JSON</a> file suitable to be
     * read by ParameterHandler::parse_input_from_json() again.
     */
    JSON = 0x0100,

    /**
     * Write the content of ParameterHandler without comments or changed default
     * values.
     */
    ShortPRM = PRM | Short,

    /**
     * Write the content of ParameterHandler without comments or changed default
     * values.
     *
     * @deprecated Use `ShortPRM` instead of `ShortText`.
     */
    ShortText = ShortPRM,

    /**
     * Write the content of ParameterHandler without comments or changed default
     * values as a XML file.
     */
    ShortXML = XML | Short,

    /**
     * Write the content of ParameterHandler without comments or changed default
     * values as a JSON file.
     */
    ShortJSON = JSON | Short,

    /**
     * Write the content of ParameterHandler without comments or changed default
     * values as a LaTeX file.
     */
    ShortLaTeX = LaTeX | Short,
  };



  /**
   * Constructor.
   */
  ParameterHandler();

  /**
   * Destructor. Declare this only to have a virtual destructor, which is
   * safer as we have virtual functions.  It actually does nothing
   * spectacular.
   */
  virtual ~ParameterHandler() override = default;

  /**
   * Inhibit automatic CopyConstructor.
   */
  ParameterHandler(const ParameterHandler &) = delete;

  /**
   * Inhibit automatic assignment operator.
   */
  ParameterHandler &
  operator=(const ParameterHandler &) = delete;

  /**
   * Parse each line from a stream until the stream returns the <tt>eof</tt>
   * condition or error to provide values for known parameter fields. The second
   * argument can be used to denote the name of the file (if that's what the
   * input stream represents) we are reading from; this is only used when
   * creating output for exceptions.
   *
   * If non-empty @p last_line is provided, the ParameterHandler object
   * will stop parsing lines after encountering @p last_line .
   * This is handy when adding extra data that shall be parsed manually.
   *
   * If @p skip_undefined is <code>true</code>, the parameter handler
   * will skip undefined sections and entries. This is useful for partially
   * parsing a parameter file, for example to obtain only the spatial dimension
   * of the problem. By default all entries and subsections are expected to be
   * declared.
   *
   * The function sets the value of all parameters it encounters in the
   * input file to the provided value. Parameters not explicitly listed
   * in the input file are left at the value they previously held, which
   * will be the default value provided to declare_entry() unless one
   * has previously read a different input file.
   *
   * Each parameter value is matched against the pattern for this
   * parameter that was provided to declare_entry(), and for each parameter
   * all associated actions that may previously have been set by
   * add_action() are executed. If a parameter does not satisfy its
   * pattern, or if an associated action throws an exception, then the
   * value provided for the parameter is not set and the current
   * object reverts to the subsection it was in before the current
   * function was called. No further processing of the input stream
   * occurs, that is everything that comes after the parameter whose
   * value does not satisfy its pattern is ignored.
   */
  virtual void
  parse_input(std::istream &     input,
              const std::string &filename       = "input file",
              const std::string &last_line      = "",
              const bool         skip_undefined = false);

  /**
   * Parse input from a specified parameter file @p filename independently
   * of the type of input file (prm, xml, json) being used. The code path
   * selected by this function is extracted from the ending of the filename,
   * so the user has to make sure that the content of the input file is
   * consistent with its name.
   *
   * The parameter @p last_line will only be used for parameter files of .prm type.
   * See the other parse_input function for documentation.
   *
   * The user can specify whether parameters in the input file not added to the
   * parameter handler will be skipped by @p skip_undefined (enables partial
   * parsing), and whether the code will assert that all parameters of the
   * parameter handler declared with flag `has_to_be_set=true` are indeed found
   * in the input file.
   *
   * If the function is called with `skip_undefined=true`, it is recommended to
   * also set `assert_mandatory_entries_are_found=true`. For example, this
   * ensures that parameters with typos in the input file will not be skipped,
   * while such mistakes would otherwise remain unrecognized.
   */
  virtual void
  parse_input(const std::string &filename,
              const std::string &last_line                          = "",
              const bool         skip_undefined                     = false,
              const bool         assert_mandatory_entries_are_found = false);

  /**
   * Parse input from a string to populate known parameter fields. The lines
   * in the string must be separated by <tt>@\n</tt> characters.
   *
   * The function in essence reads the entire file into a stream and
   * then calls the other parse_input() function with that stream. See
   * there for more information.
   */
  virtual void
  parse_input_from_string(const std::string &s,
                          const std::string &last_line      = "",
                          const bool         skip_undefined = false);

  /**
   * Parse input from an XML stream to populate known parameter fields. This
   * could be from a file originally written by the print_parameters() function
   * using the XML output style and then modified by hand as necessary, or from
   * a file written using this method and then modified by the graphical
   * parameter GUI (see the general documentation of this class).
   */
  virtual void
  parse_input_from_xml(std::istream &input, const bool skip_undefined = false);

  /**
   * Parse input from a JSON stream to populate known parameter fields. This
   * could be from a file originally written by the print_parameters() function
   * using the JSON output style and then modified by hand as necessary, or from
   * a separate program that knows how to write JSON format for ParameterHandler
   * input.
   */
  virtual void
  parse_input_from_json(std::istream &input, const bool skip_undefined = false);

  /**
   * Clear all contents.
   */
  void
  clear();


  /**
   * Declare a new entry with name <tt>entry</tt>, default and for which any
   * input has to match the <tt>pattern</tt> (default: any pattern).
   *
   * The function generates an exception of type ExcValueDoesNotMatchPattern
   * if the default value doesn't match the given pattern, using the C++ throw
   * mechanism. However, this exception is only generated <i>after</i> the
   * entry has been created; if you have code where no sensible default value
   * for a parameter is possible, you can then catch and ignore this
   * exception.
   *
   * The parameter @p documentation defaulting to an empty string is used
   * to add a documenting text to each entry which will be printed as a comment
   * when this class is asked to write out all declarations to a stream using
   * the print_parameters() function.
   *
   * The parameter @p has_to_be_set can be used in order to declare this
   * parameter as a parameter whose default value has to be overwritten by
   * one of the methods provided by this class. Whether a parameter has been set
   * successfully can be queried by the functions get_entries_wrongly_not_set()
   * and assert_that_entries_have_been_set().
   *
   * @note An entry can be declared more than once without generating an
   * error, for example to override an earlier default value.
   */
  void
  declare_entry(const std::string &          entry,
                const std::string &          default_value,
                const Patterns::PatternBase &pattern = Patterns::Anything(),
                const std::string &          documentation = "",
                const bool                   has_to_be_set = false);

  /**
   * Attach an action to the parameter with name @p entry in the current
   * section. The action needs to be a function-like object that takes the
   * value of the parameter as a (string) argument. See the general
   * documentation of this class for a longer description of actions, as well as
   * examples.
   *
   * The action is executed in three different circumstances:
   * - With the default value of the parameter with name @p name, at
   *   the end of the current function. This is useful because it allows
   *   for the action to execute whatever it needs to do at least once
   *   for each parameter, even those that are not actually specified in
   *   the input file (and thus remain at their default values).
   * - Within the ParameterHandler::set() functions that explicitly
   *   set a value for a parameter.
   * - Within the parse_input() function and similar functions such
   *   as parse_input_from_string(). Here, the action is executed
   *   whenever the parameter with which it is associated is read
   *   from the input, after it has been established that the value
   *   so read matches the pattern that corresponds to this parameter,
   *   and before the value is actually saved.
   *
   * It is valid to add multiple actions to the same parameter. They will
   * in that case be executed in the same order in which they were added.
   *
   * @note Actions may modify all sorts of variables in their scope. The
   *  only thing an action should not modify is the ParameterHandler object
   *  it is attached to. In other words, it is not allowed to enter or
   *  leave sections of the current ParameterHandler object. It is, in
   *  principle, acceptable to call ParameterHandler::get() and related
   *  functions on other parameters in the current section, but since
   *  there is no guarantee about the order in which they will be read
   *  from an input file, you will not want to rely on the values these
   *  functions would return.
   *
   * @note Throwing an exception in an action is generally not a good
   *  idea, but yields fundamentally the same result as if one tries to
   *  read a parameter from a file for which the value does not satisfy
   *  the pattern associated with the parameter. In other words, the
   *  value just read is discarded, and ParameterHandler::parse_input()
   *  stops to read any further content from the file. See
   *  ParameterHandler::parse_input() for more information.
   */
  void
  add_action(const std::string &                                  entry,
             const std::function<void(const std::string &value)> &action);

  /**
   * Declare a new entry name @p entry, set its default value to the content of
   * the variable @p parameter, and create an action that will fill @p
   * parameter with updated values when a file is parsed, or the entry is set
   * to a new value.
   *
   * By default, the pattern to use is obtained by calling the function
   * Patterns::Tools::Convert<T>::to_pattern(), but a custom one can be used.
   *
   * The parameter @p has_to_be_set can be used in order to declare this
   * parameter as a parameter whose default value has to be overwritten by
   * one of the methods provided by this class. Whether a parameter has been set
   * successfully can be queried by the functions get_entries_wrongly_not_set()
   * and assert_that_entries_have_been_set().
   */
  template <class ParameterType>
  void
  add_parameter(const std::string &          entry,
                ParameterType &              parameter,
                const std::string &          documentation = "",
                const Patterns::PatternBase &pattern =
                  *Patterns::Tools::Convert<ParameterType>::to_pattern(),
                const bool has_to_be_set = false);

  /**
   * Create an alias for an existing entry. This provides a way to refer to a
   * parameter in the input file using an alternate name. The alias will be in
   * the current section, and the referenced entry needs to be an existing
   * entry in the current section.
   *
   * The primary purpose of this function is to allow for a backward
   * compatible way of changing names in input files of applications for which
   * backward compatibility is important. This can be achieved by changing the
   * name of the parameter in the call to declare_entry(), and then creating
   * an alias that maps the old name to the new name. This way, old input
   * files can continue to refer to parameters under the old name, and they
   * will automatically be mapped to the new parameter name.
   *
   * It is valid to set the same parameter multiple times in an input file.
   * The value that will ultimately be chosen in such cases is simply the last
   * value set. This rule also applies to aliases, where the final value of a
   * parameter is the last value set either through the current name of the
   * parameter or through any of its possible multiple aliases. For example,
   * if you have an input file that looks like
   * @code
   *   set parm1       = 1
   *   set parm1_alias = 2
   * @endcode
   * where <code>parm1_alias</code> is an alias declared via
   * @code
   *   prm.declare_alias ("parm1", "parm1_alias");
   * @endcode
   * then the final value for the parameter called <code>parm1</code> will be
   * 2, not 1.
   *
   * @param existing_entry_name The name of an existing parameter in the
   * current section that the alias should refer to.
   * @param alias_name An alternate name for the parameter referenced by the
   * first argument.
   * @param alias_is_deprecated If true, mark the alias as deprecated. This
   * will then be listed in the description of the alias if you call
   * print_parameters(), and you will get a warning on the screen when reading
   * an input file that contains this deprecated alias. The purpose of this
   * argument is to be able to allow the use of an old name for a parameter
   * (see above) but make it clear that this old name will eventually be
   * removed.
   */
  void
  declare_alias(const std::string &existing_entry_name,
                const std::string &alias_name,
                const bool         alias_is_deprecated = false);

  /**
   * Enter a subsection. If it does not yet exist, create it.
   */
  void
  enter_subsection(const std::string &subsection);

  /**
   * Leave present subsection.
   */
  void
  leave_subsection();

  /**
   * Check whether a subsection or a subsection path exists in current tree.
   * The input parameter @p sub_path is assumed to be relative to the
   * currently selected path.
   */
  bool
  subsection_path_exists(const std::vector<std::string> &sub_path) const;

  /**
   * Return value of entry @p entry_string.  If the entry was changed,
   * then the changed value is returned, otherwise the default value. If the
   * value of an undeclared entry is required, an @p Assert will fail.
   */
  std::string
  get(const std::string &entry_string) const;

  /**
   * Return value of entry @p entry_string.  If the entry was changed,
   * then the changed value is returned, otherwise the default value. If the
   * value of an undeclared entry is required, an @p Assert will fail.
   * If @p entry_subsection_path is non-empty, the value will be gotten
   * from the subsection represented by that path instead of the current
   * subsection. The first string in @p entry_subsection_path must be the name
   * of a subsection of the current section, and each next string must be the
   * name of a subsection of the one before it.
   */
  std::string
  get(const std::vector<std::string> &entry_subsection_path,
      const std::string &             entry_string) const;

  /**
   * Return value of entry @p entry_string as <code>long int</code>. (A long
   * int is chosen so that even very large unsigned values can be returned by
   * this function).
   */
  long int
  get_integer(const std::string &entry_string) const;

  /**
   * Return value of entry @p entry_string as <code>long int</code>. (A long
   * int is chosen so that even very large unsigned values can be returned by
   * this function).
   * If @p entry_subsection_path is non-empty, the value will be gotten
   * from the subsection represented by that path instead of the current
   * subsection.
   */
  long int
  get_integer(const std::vector<std::string> &entry_subsection_path,
              const std::string &             entry_string) const;

  /**
   * Return value of entry @p entry_name as @p double.
   */
  double
  get_double(const std::string &entry_name) const;

  /**
   * Return value of entry @p entry_name as @p double.
   * If @p entry_subsection_path is non-empty, the value will be gotten
   * from the subsection represented by that path instead of the current
   * subsection.
   */
  double
  get_double(const std::vector<std::string> &entry_subsection_path,
             const std::string &             entry_string) const;
  /**
   * Return value of entry @p entry_name as @p bool. The entry may
   * be "true" or "yes" for @p true, "false" or "no" for @p false
   * respectively.
   */
  bool
  get_bool(const std::string &entry_name) const;

  /**
   * Return value of entry @p entry_name as @p bool. The entry may
   * be "true" or "yes" for @p true, "false" or "no" for @p false
   * respectively.
   * If @p entry_subsection_path is non-empty, the value will be gotten
   * from the subsection represented by that path instead of the current
   * subsection.
   */
  bool
  get_bool(const std::vector<std::string> &entry_subsection_path,
           const std::string &             entry_string) const;

  /**
   * Change the value presently stored for <tt>entry_name</tt> to the one
   * given in the second argument.
   *
   * The parameter must already exist in the present subsection.
   *
   * The function throws an exception of type ExcValueDoesNotMatchPattern if
   * the new value does not conform to the pattern for this entry.
   */
  void
  set(const std::string &entry_name, const std::string &new_value);

  /**
   * Same as above, but an overload where the second argument is a character
   * pointer. This is necessary, since otherwise the call to
   * <tt>set("abc","def")</code> will be mapped to the function taking one
   * string and a bool as arguments, which is certainly not what is most often
   * intended.
   *
   * The function throws an exception of type ExcValueDoesNotMatchPattern if
   * the new value does not conform to the pattern for this entry.
   */
  void
  set(const std::string &entry_name, const char *new_value);

  /**
   * Change the value presently stored for <tt>entry_name</tt> to the one
   * given in the second argument.
   *
   * The parameter must already exist in the present subsection.
   *
   * The function throws an exception of type ExcValueDoesNotMatchPattern if
   * the new value does not conform to the pattern for this entry.
   */
  void
  set(const std::string &entry_name, const long int new_value);

  /**
   * Change the value presently stored for <tt>entry_name</tt> to the one
   * given in the second argument.
   *
   * The parameter must already exist in the present subsection.
   *
   * For internal purposes, the new value needs to be converted to a string.
   * This is done using 16 digits of accuracy, so the set value and the one
   * you can get back out using get_double() may differ in the 16th digit.
   *
   * The function throws an exception of type ExcValueDoesNotMatchPattern if
   * the new value does not conform to the pattern for this entry.
   */
  void
  set(const std::string &entry_name, const double new_value);

  /**
   * Change the value presently stored for <tt>entry_name</tt> to the one
   * given in the second argument.
   *
   * The parameter must already exist in the present subsection.
   *
   * The function throws an exception of type ExcValueDoesNotMatchPattern if
   * the new value does not conform to the pattern for this entry.
   */
  void
  set(const std::string &entry_name, const bool new_value);

  /**
   * Print all parameters with the given @p style to @p out.
   *
   * Before printing, all current parameters and subsections are sorted
   * alphabetically by default.
   * This behavior can be disabled setting the optional parameter @p style
   * to <tt>KeepDeclarationOrder</tt>: in this case entries are printed in the
   * same order as they have been declared.
   *
   * In <tt>PRM</tt>, <tt>XML</tt>, and <tt>JSON</tt> format, the output is
   * formatted in such a way that it is
   * possible to use it for later input again. This is most useful to record
   * the parameters for a specific run, since if you output the parameters
   * using this function into a log file, you can always recover the results
   * by simply copying the output to your input file.
   *
   * Besides the name and value of each entry, the output also contains the
   * default value of entries if it is different from the actual value, as
   * well as the documenting string given to the declare_entry() function if
   * available.
   *
   * By using the flag <tt>Short</tt> in combination with <tt>PRM</tt>,
   * <tt>XML</tt>, <tt>JSON</tt>, or <tt>LaTeX</tt> (or by using the shortcuts
   * <tt>ShortPRM</tt>, <tt>ShortXML</tt>, <tt>ShortJSON</tt>, or
   * <tt>ShortLaTeX</tt>), a reduced output can be generated, only containing
   * the values and skipping the documentation.
   *
   * In <tt>XML</tt> format, the output starts with one root element
   * <tt>ParameterHandler</tt> in order to get a valid XML document and all
   * subsections under it.
   *
   * In <tt>LaTeX</tt> format, the output contains the same information but in
   * a format so that the resulting file can be input into a latex document
   * such as a manual for the code for which this object handles run-time
   * parameters. The various sections of parameters are then represented by
   * latex section and subsection commands as well as by nested enumerations.
   *
   * You can reference specific parameter sections and individual parameters
   * by the labels that are generated automatically for each entry. The
   * labels have the format <code>parameters:section1/subsection1</code>
   * and <code>parameters:section1/subsection1/someentry</code>. Because
   * special characters can appear in the section and entry names, these
   * will be "mangled". Here, all characters except <code>[a-zA-Z0-9]</code>
   * are replaced by <code>_XX</code>, where <code>XX</code> is the two-digit
   * ascii code of the character in hexadecimal encoding (so a space becomes
   * <code>_20</code> for example).
   *
   * While this function escapes special LaTeX-specific characters (backslash,
   * underscore, etc.) in most of the output (names, default values, etc.),
   * the documentation string is passed as-is. This means you can use math
   * environments and other formatting in the description, but you need
   * to escape quotes, backslashes, underscores, etc. yourself.
   *
   * In addition, all parameter names are listed with <code>@\index</code>
   * statements in two indices called <code>prmindex</code> (where the name of
   * each parameter is listed in the index) and <code>prmindexfull</code>
   * where parameter names are listed sorted by the section in which they
   * exist. By default, the LaTeX program ignores these <code>@\index</code>
   * commands, but they can be used to generate an index by using the
   * following commands in the preamble of the latex file:
   * @code
   * \usepackage{imakeidx}
   * \makeindex[name=prmindex, title=Index of run-time parameter entries]
   * \makeindex[name=prmindexfull,
   *            title=Index of run-time parameters with section names]
   * @endcode
   * and at the end of the file this:
   * @code
   * \printindex[prmindex]
   * \printindex[prmindexfull]
   * @endcode
   */
  std::ostream &
  print_parameters(std::ostream &out, const OutputStyle style) const;



  /**
   * Print all parameters to the file given by @p filename with the given output
   * style @p style.
   *
   * This function deduces the output format from the extension of the specified
   * filename. Supported extensions are `prm`, `xml`, `tex`, and `json`. Hence,
   * it is not necessary to specify an output format via the @p style argument
   * as long as one of these extensions is added to the filename. If an output
   * format is specified in the @p style parameter nevertheless, the output
   * format has to be consistent with the filename extension.
   *
   * If no extension is specified or the extension is not supported, the
   * output format is deduced from the @p style argument.
   *
   * If neither the extension is supported, nor does the @p style parameter
   * contain a format specification, an assertion is thrown.
   *
   * @param filename The output file name.
   * @param style The style with which output is produced.
   */
  void
  print_parameters(const std::string &filename,
                   const OutputStyle  style = DefaultStyle) const;

  /**
   * Print parameters to a logstream. This function allows to print all
   * parameters into a log-file. Sections will be indented in the usual log-
   * file style.
   *
   * All current parameters and subsections are sorted
   * alphabetically by default.
   * This behavior can be disabled setting the optional parameter @p style
   * to <tt>KeepDeclarationOrder</tt>: in this case entries are printed in the
   * same order as they have been declared.
   *
   * @note All style settings in @p style not related to the ordering are
   *   ignored.
   */
  void
  log_parameters(LogStream &out, const OutputStyle style = DefaultStyle);

  /**
   * Log parameters in the present subsection. The subsection is determined by
   * the <tt>subsection_path</tt> member variable. This variable is controlled
   * by entering and leaving subsections through the enter_subsection() and
   * leave_subsection() functions.
   *
   * All current parameters and subsections are sorted
   * alphabetically by default.
   * This behavior can be disabled setting the optional parameter @p style
   * to <tt>KeepDeclarationOrder</tt>: in this case entries are printed in the
   * same order as they have been declared.
   *
   * @note All style settings in @p style not related to the ordering are
   *   ignored.
   *
   * In most cases, you will not want to use this function directly, but have
   * it called recursively by the previous function.
   */
  void
  log_parameters_section(LogStream &       out,
                         const OutputStyle style = DefaultStyle);

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Write the data of this object to a stream for the purpose of
   * serialization.
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * Read the data of this object from a stream for the purpose of
   * serialization.
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * Write and read the data of this object from a stream for the purpose
   * of serialization.
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  /**
   * Test for equality.
   */
  bool
  operator==(const ParameterHandler &prm2) const;

  /**
   * Return a set of parameter names (including subsection names) corresponding
   * to those entries of the parameter handler that have not been set by one of
   * the functions parsing parameters from an input file or by an explicit call
   * to one of the set() functions, but that have been declared as mandatory
   * parameters that must be set (through the last argument of the
   * declare_entry() function or add_parameter() function).
   */
  std::set<std::string>
  get_entries_wrongly_not_set() const;

  /**
   * Asserts that those entries of the parameter handler with flag
   * `has_to_be_set = true` have been set. An exception is invoked
   * if at least one of these parameters has not been set.
   */
  void
  assert_that_entries_have_been_set() const;

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception
   */
  DeclException1(ExcEntryAlreadyExists,
                 std::string,
                 << "The following entry already exists: " << arg1 << ".");
  /**
   * Exception
   */
  DeclException2(ExcValueDoesNotMatchPattern,
                 std::string,
                 std::string,
                 << "The string <" << arg1
                 << "> does not match the given pattern <" << arg2 << ">.");
  /**
   * Exception
   */
  DeclExceptionMsg(
    ExcAlreadyAtTopLevel,
    "You can't leave a subsection if you are already at the top level "
    "of the subsection hierarchy.");
  /**
   * Exception
   */
  DeclException1(ExcEntryUndeclared,
                 std::string,
                 << "You can't ask for entry <" << arg1
                 << "> you have not yet declared.");

  /**
   * Exception for when there are an unequal number of 'subsection' and 'end'
   * statements. The first argument is the name of the file and the second
   * argument is a formatted list of the subsection path before and after
   * entering the parser.
   */
  DeclException2(ExcUnbalancedSubsections,
                 std::string,
                 std::string,
                 << "There are unequal numbers of 'subsection' and 'end' "
                    "statements in the parameter file <"
                 << arg1 << ">." << (arg2.size() > 0 ? "\n" + arg2 : ""));

  /**
   * Exception for when, during parsing of a parameter file, the parser
   * encounters a subsection in the file that was not previously declared.
   */
  DeclException3(ExcNoSubsection,
                 int,
                 std::string,
                 std::string,
                 << "Line <" << arg1 << "> of file <" << arg2
                 << ": There is "
                    "no such subsection to be entered: "
                 << arg3);

  /**
   * General exception for a line that could not be parsed, taking, as
   * arguments, the line number, file name, and a brief description of why the
   * line cannot be parsed.
   */
  DeclException3(ExcCannotParseLine,
                 int,
                 std::string,
                 std::string,
                 << "Line <" << arg1 << "> of file <" << arg2 << ">: " << arg3);

  /**
   * Exception for an entry in a parameter file that does not match the
   * provided pattern. The arguments are, in order, the line number, file
   * name, entry value, entry name, and a description of the pattern.
   */
  DeclException5(ExcInvalidEntryForPattern,
                 int,
                 std::string,
                 std::string,
                 std::string,
                 std::string,
                 << "Line <" << arg1 << "> of file <" << arg2
                 << ">:\n"
                    "    The entry value \n"
                 << "        " << arg3 << '\n'
                 << "    for the entry named\n"
                 << "        " << arg4 << '\n'
                 << "    does not match the given pattern:\n"
                 << "        " << arg5);

  /**
   * Exception for when an XML file cannot be read at all. This happens when
   * there is no top-level XML element called "ParameterHandler" or when there
   * are multiple top level elements.
   */
  DeclExceptionMsg(ExcInvalidXMLParameterFile,
                   "The provided file could not be parsed as a "
                   "ParameterHandler description.");

  /**
   * Exception for when an entry in an XML parameter file does not match the
   * provided pattern. The arguments are, in order, the entry value, entry
   * name, and a description of the pattern.
   *
   * @deprecated Use ExcValueDoesNotMatchPattern instead of ExcInvalidEntryForPatternXML.
   */
  using ExcInvalidEntryForPatternXML DEAL_II_DEPRECATED =
    ExcValueDoesNotMatchPattern;

  /**
   * Exception for when the file given in an include statement cannot be
   * open. The arguments are, in order, the line number of the include
   * statement, current parameter file name, and the name of the file intended
   * for inclusion.
   */
  DeclException3(
    ExcCannotOpenIncludeStatementFile,
    int,
    std::string,
    std::string,
    << "Line <" << arg1 << "> of file <" << arg2
    << ">: This line "
       "contains an 'include' or 'INCLUDE' statement, but the given "
       "file to include <"
    << arg3 << "> cannot be opened.");

  //@}

private:
  /**
   * The separator used when accessing elements of a path into the parameter
   * tree.
   */
  static const char path_separator = '.';

  /**
   * Path of presently selected subsections; empty list means top level
   */
  std::vector<std::string> subsection_path;

  /**
   * The complete tree of sections and entries. See the general documentation
   * of this class for a description how data is stored in this variable.
   *
   * The variable is a pointer so that we can use an incomplete type, rather
   * than having to include all of the property_tree stuff from boost. This
   * works around a problem with gcc 4.5.
   */
  std::unique_ptr<boost::property_tree::ptree> entries;

  /**
   * A map that stores a pair of boolean variables for each entry
   * added to the parameter handler. The first bool describes whether the
   * parameter has to be set according to the last argument of the functions
   * declare_entry() or add_parameter(), and the second bool contains the
   * information whether the parameter has been set by any of the functions
   * parsing input parameters or by a set function of this class.
   */
  std::map<std::string, std::pair<bool, bool>> entries_set_status;

  /**
   * A list of patterns that are used to describe the parameters of this
   * object. Every nodes in the property tree corresponding to a parameter
   * stores an index into this array.
   */
  std::vector<std::unique_ptr<const Patterns::PatternBase>> patterns;

  /**
   * A list of actions that are associated with parameters. These
   * are added by the add_action() function. Nodes in the property
   * tree corresponding to individual parameters
   * store indices into this array in order to reference specific actions.
   */
  std::vector<std::function<void(const std::string &)>> actions;

  /**
   * Return the string that identifies the current path into the property
   * tree. This is only a path, i.e. it is not terminated by the
   * path_separator character.
   *
   * This function simply calls collate_path_string() with
   * @p subsection_path as argument
   */
  std::string
  get_current_path() const;

  /**
   * Given the name of an entry as argument, the function computes a full path
   * into the parameter tree using the current subsection.
   */
  std::string
  get_current_full_path(const std::string &name) const;

  /**
   * This function computes a full path into the parameter tree given a path
   * from the current subsection and the name of an entry.
   */
  std::string
  get_current_full_path(const std::vector<std::string> &sub_path,
                        const std::string &             name) const;

  /**
   * Scan one line of input. <tt>input_filename</tt> and
   * <tt>current_line_n</tt> are the name of the input file and the number of
   * the line presently scanned (these are used in exception messages to show
   * where parse errors occurred). This function will raise an exception if
   * the line contains an undeclared subsection or entry, if the line's entry
   * does not match its given pattern, or if the line could not be understood
   * as a valid parameter file expression.
   *
   * The function modifies its argument, but also takes it by value, so the
   * caller's variable is not changed.
   *
   * If @p skip_undefined is <code>true</code>, the parser
   * will skip undefined sections and entries. This is useful for partially
   * parsing a parameter file, for example to obtain only the spatial dimension
   * of the problem. By default all entries and subsections are expected to be
   * declared.
   */
  void
  scan_line(std::string        line,
            const std::string &input_filename,
            const unsigned int current_line_n,
            const bool         skip_undefined);

  /**
   * Print out the parameters of the subsection given by the
   * @p target_subsection_path argument, as well as all subsections
   * within it recursively. This function is called from the
   * print_parameters() function, and is implemented for all @p style
   * arguments other than XML and JSON (where we can output the
   * entire set of parameters via BOOST functions). The @p indent_level
   * argument indicates how many spaces the output should be indented,
   * so that subsections properly nest inside the output of higher
   * sections.
   */
  void
  recursively_print_parameters(
    const boost::property_tree::ptree & tree,
    const std::vector<std::string> &    target_subsection_path,
    const ParameterHandler::OutputStyle style,
    const unsigned int                  indent_level,
    std::ostream &                      out) const;

  friend class MultipleParameterLoop;
};

/**
 * Global operator which returns an object in which all bits are set which are
 * either set in the first or the second argument. This operator exists since
 * if it did not then the result of the bit-or <tt>operator |</tt> would be an
 * integer which would in turn trigger a compiler warning when we tried to
 * assign it to an object of type ParameterHandler::OutputStyle.
 */
inline ParameterHandler::OutputStyle
operator|(const ParameterHandler::OutputStyle f1,
          const ParameterHandler::OutputStyle f2)
{
  return static_cast<ParameterHandler::OutputStyle>(
    static_cast<unsigned int>(f1) | static_cast<unsigned int>(f2));
}

/**
 * The class MultipleParameterLoop offers an easy possibility to test several
 * parameter sets during one run of the program. For this it uses the
 * ParameterHandler class to read in data in a standardized form, searches for
 * variant entry values and performs a loop over all combinations of
 * parameters.
 *
 * Variant entry values are given like this:
 *   @verbatim
 *     set Time step size = { 0.1 | 0.2 | 0.3 }
 *   @endverbatim
 * The loop will then perform three runs of the program, one for each value of
 * <tt>Time step size</tt>, while all other parameters are as specified or
 * with their default value. If there are several variant entry values in the
 * input, a loop is performed for each combination of variant values:
 *   @verbatim
 *     set Time step size = { 0.1 | 0.2 }
 *     set Solver         = { CG  | GMRES }
 *   @endverbatim
 * will result in four runs of the programs, with time step 0.1 and 0.2 for
 * each of the two solvers.
 *
 * In addition to variant entries, this class also supports <i>array
 * entries</i> that look like this:
 *   @verbatim
 *     set Output file = ofile.{{ 1 | 2 | 3 | 4 }}
 *   @endverbatim
 * This indicates that if there are variant entries producing a total of four
 * different runs, then we will write their results to the files
 * <tt>ofile.1</tt>, <tt>ofile.2</tt>, <tt>ofile.3</tt> and <tt>ofile.4</tt>,
 * respectively. Array entries do not generate multiple runs of the main loop
 * themselves, but if there are variant entries, then in the <i>n</i>th run of
 * the main loop, also the <i>n</i>th value of an array is returned.
 *
 * Since the different variants are constructed in the order of declaration,
 * not in the order in which the variant entries appear in the input file, it
 * may be difficult to guess the mapping between the different variants and
 * the appropriate entry in an array. You will have to check the order of
 * declaration, or use only one variant entry.
 *
 * It is guaranteed that only selections which match the regular expression
 * (pattern) given upon declaration of an entry are given back to the program.
 * If a variant value does not match the regular expression, the default value
 * is stored and an error is issued. Before the first run of the loop, all
 * possible values are checked for their conformance, so that the error is
 * issued at the very beginning of the program.
 *
 *
 * <h3>Usage</h3>
 *
 * The usage of this class is similar to the ParameterHandler class. First the
 * entries and subsections have to be declared, then a loop is performed in
 * which the different parameter sets are set, a new instance of a user class
 * is created which is then called. Taking the classes of the example for the
 * ParameterHandler class, the extended program would look like this:
 *   @code
 *     class HelperClass : public MultipleParameterLoop::UserClass
 *     {
 *     public:
 *       HelperClass ();
 *
 *       virtual void create_new (const unsigned int run_no);
 *       virtual void run (ParameterHandler &prm);
 *
 *       static void declare_parameters (ParameterHandler &prm);
 *     private:
 *       std::unique_ptr<Problem> p;
 *     };
 *
 *
 *     HelperClass::HelperClass () : p(0) {}
 *
 *
 *     void HelperClass::create_new (const unsigned int run_no)
 *     {
 *       p = std::make_unique<Problem>());
 *     }
 *
 *
 *     void HelperClass::declare_parameters (ParameterHandler &prm)
 *     {
 *       Problem::declare_parameters (prm);
 *     }
 *
 *
 *     void HelperClass::run (ParameterHandler &prm)
 *     {
 *       p->get_parameters (prm);
 *       p->do_useful_work ();
 *     }
 *
 *
 *
 *     int main ()
 *     {
 *       class MultipleParameterLoop prm;
 *       HelperClass h;
 *       HelperClass::declare_parameters (prm);
 *       prm.parse_input ("prmtest.prm");
 *       prm.loop (h);
 *       return 0;
 *     }
 *   @endcode
 *
 * As can be seen, first a new helper class has to be set up. This must
 * contain a virtual constructor for a problem class. You can also derive your
 * problem class from MultipleParameterLoop::UserClass and let
 * <tt>create_new</tt> clear all member variables. If you have access to all
 * inherited member variables in some way this is the recommended procedure. A
 * third possibility is to use multiple inheritance and derive a helper class
 * from both the MultipleParameterLoop::UserClass and the problem class. In
 * any case, <tt>create_new</tt> has to provide a clean problem object which
 * is the problem in the second and third possibility.
 *
 * The derived class also has to provide for member functions which declare
 * the entries and which run the program. Running the program includes getting
 * the parameters out of the ParameterHandler object.
 *
 * After defining an object of this helper class and an object of the
 * MultipleParameterLoop class, the entries have to be declared in the same
 * way as for the ParameterHandler class. Then the input has to be read.
 * Finally the loop is called. This executes the following steps:
 *   @code
 *     for (each combination)
 *       {
 *         UserObject.create_new (run_no);
 *
 *         // set parameters for this run
 *
 *         UserObject.run (*this);
 *       }
 *   @endcode
 * <tt>UserObject</tt> is the parameter to the <tt>loop</tt> function.
 * <tt>create_new</tt> is given the number of the run (starting from one) to
 * enable naming output files differently for each run.
 *
 *
 * <h3>Syntax for variant and array entry values</h3>
 *
 * Variant values are specified like <tt>prefix{ v1 | v2 | v3 | ...
 * }postfix</tt>. Whitespace to the right of the opening brace <tt>{</tt> is
 * ignored as well as to the left of the closing brace <tt>}</tt> while
 * whitespace on the respectively other side is not ignored. Whitespace around
 * the mid symbols <tt>|</tt> is also ignored. The empty selection <tt>prefix{
 * v1 | }postfix</tt> is also allowed and produces the strings
 * <tt>prefixv1postfix</tt> and <tt>prefixpostfix</tt>.
 *
 * The syntax for array values is equal, apart from the double braces:
 * <tt>prefix{{ v1 | v2 | v3 }}postfix</tt>.
 *
 *
 * <h3>Worked example</h3>
 *
 * Given the above extensions to the example program for the ParameterHandler
 * and the following input file
 *   @verbatim
 *     set Equation 1 = Poisson
 *     set Equation 2 = Navier-Stokes
 *     set Output file= results.{{ 1 | 2 | 3 | 4 | 5 | 6 }}
 *
 *     subsection Equation 1
 *       set Matrix type = Sparse
 *       subsection Linear solver
 *         set Solver                       = CG
 *         set Maximum number of iterations = { 10 | 20 | 30 }
 *       end
 *     end
 *
 *     subsection Equation 2
 *       set Matrix type = Full
 *       subsection Linear solver
 *         set Solver                       = { BiCGStab | GMRES }
 *         set Maximum number of iterations = 100
 *       end
 *     end
 *   @endverbatim
 * this is the output:
 *   @verbatim
 *     LinEq: Method=CG, MaxIterations=10
 *     LinEq: Method=BiCGStab, MaxIterations=100
 *     Problem: outfile=results.1
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=20
 *     LinEq: Method=BiCGStab, MaxIterations=100
 *     Problem: outfile=results.2
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=30
 *     LinEq: Method=BiCGStab, MaxIterations=100
 *     Problem: outfile=results.3
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=10
 *     LinEq: Method=GMRES, MaxIterations=100
 *     Problem: outfile=results.4
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=20
 *     LinEq: Method=GMRES, MaxIterations=100
 *     Problem: outfile=results.5
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *     LinEq: Method=CG, MaxIterations=30
 *     LinEq: Method=GMRES, MaxIterations=100
 *     Problem: outfile=results.6
 *              eq1=Poisson, eq2=Navier-Stokes
 *              Matrix1=Sparse, Matrix2=Full
 *   @endverbatim
 * Since <tt>create_new</tt> gets the number of the run it would also be
 * possible to output the number of the run.
 *
 *
 * @ingroup input
 * @author Wolfgang Bangerth, October 1997, 2010
 */
class MultipleParameterLoop : public ParameterHandler
{
public:
  /**
   * This is the class the helper class or the problem class has to be derived
   * of.
   */
  class UserClass
  {
  public:
    /**
     * Destructor. It doesn't actually do anything, but is declared to force
     * derived classes to have a virtual destructor.
     */
    virtual ~UserClass() = default;

    /**
     * <tt>create_new</tt> must provide a clean object, either by creating a
     * new one or by cleaning an old one.
     */
    virtual void
    create_new(const unsigned int run_no) = 0;

    /**
     * Get the parameters and run any necessary action.
     */
    virtual void
    run(ParameterHandler &prm) = 0;
  };

  /**
   * Constructor
   */
  MultipleParameterLoop();

  /**
   * Destructor. Declare this only to have a virtual destructor, which is
   * safer as we have virtual functions. It actually does nothing spectacular.
   */
  virtual ~MultipleParameterLoop() override = default;

  /**
   * Read input from a stream until the stream returns the <tt>eof</tt>
   * condition or error. The second argument can be used to denote the name of
   * the file (if that's what the input stream represents) we are reading
   * from; this is only used when creating output for error messages.
   *
   * If non-empty @p last_line is provided, the ParameterHandler object
   * will stop parsing lines after encountering @p last_line .
   * This is handy when adding extra data that shall be parsed manually.
   *
   * If @p skip_undefined is <code>true</code>, the parameter handler
   * will skip undefined sections and entries. This is useful for partially
   * parsing a parameter file, for example to obtain only the spatial dimension
   * of the problem. By default all entries and subsections are expected to be
   * declared.
   *
   * @note This is the only overload of the three <tt>parse_input</tt>
   * functions implemented by ParameterHandler overridden with new behavior by
   * this class. This is because the other two <tt>parse_input</tt> functions
   * just reformat their inputs and then call this version.
   */
  virtual void
  parse_input(std::istream &     input,
              const std::string &filename       = "input file",
              const std::string &last_line      = "",
              const bool         skip_undefined = false) override;

  /**
   * Overriding virtual functions which are overloaded (like
   * ParameterHandler::parse_input, which has two different sets of input
   * argument types) causes the non-overridden functions to be hidden. Get
   * around this by explicitly using both variants of
   * ParameterHandler::parse_input and then overriding the one we care about.
   */
  using ParameterHandler::parse_input;

  /**
   * run the central loop.
   */
  void
  loop(UserClass &uc);

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * An object in the list of entries with multiple values.
   */
  class Entry
  {
  public:
    /**
     * Declare what a multiple entry is: a variant entry (in curly braces
     * <tt>{</tt>, <tt>}</tt>) or an array (in double curly braces
     * <tt>{{</tt>, <tt>}}</tt>).
     */
    enum MultipleEntryType
    {
      /**
       * A variant entry.
       */
      variant,
      /**
       * An array entry.
       */
      array
    };

    /**
     * Constructor
     */
    Entry()
      : type(array)
    {}

    /**
     * Construct an object with given subsection path, name and value. The
     * splitting up into the different variants is done later by
     * <tt>split_different_values</tt>.
     */
    Entry(const std::vector<std::string> &Path,
          const std::string &             Name,
          const std::string &             Value);

    /**
     * Split the entry value into the different branches.
     */
    void
    split_different_values();

    /**
     * Path to variant entry.
     */
    std::vector<std::string> subsection_path;

    /**
     * Name of entry.
     */
    std::string entry_name;

    /**
     * Original variant value.
     */
    std::string entry_value;

    /**
     * List of entry values constructed out of what was given in the input
     * file.
     */
    std::vector<std::string> different_values;

    /**
     * Store whether this entry is a variant entry or an array.
     */
    MultipleEntryType type;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;
  };

  /**
   * List of variant entry values.
   */
  std::vector<Entry> multiple_choices;

  /**
   * Number of branches constructed from the different combinations of the
   * variants. This obviously equals the number of runs to be performed.
   */
  unsigned int n_branches;

  /**
   * Initialize the different branches, i.e.  construct the combinations.
   */
  void
  init_branches();

  /**
   * Traverse the section currently set by
   * enter_subsection()/leave_subsection() and see which of the entries are
   * variant or array entries. Then fill the multiple_choices variable using
   * this information.
   */
  void
  init_branches_current_section();

  /**
   * Transfer the entry values for one run to the entry tree.
   */
  void
  fill_entry_values(const unsigned int run_no);
};


// ---------------------- inline and template functions --------------------
template <class Archive>
inline void
ParameterHandler::save(Archive &ar, const unsigned int) const
{
  // Forward to serialization
  // function in the base class.
  ar &static_cast<const Subscriptor &>(*this);

  ar &*entries.get();

  std::vector<std::string> descriptions;

  for (const auto &pattern : patterns)
    descriptions.push_back(pattern->description());

  ar &descriptions;
}


template <class Archive>
inline void
ParameterHandler::load(Archive &ar, const unsigned int)
{
  // Forward to serialization
  // function in the base class.
  ar &static_cast<Subscriptor &>(*this);

  ar &*entries.get();

  std::vector<std::string> descriptions;
  ar &                     descriptions;

  patterns.clear();
  for (const auto &description : descriptions)
    patterns.push_back(Patterns::pattern_factory(description));
}


template <class ParameterType>
void
ParameterHandler::add_parameter(const std::string &          entry,
                                ParameterType &              parameter,
                                const std::string &          documentation,
                                const Patterns::PatternBase &pattern,
                                const bool                   has_to_be_set)
{
  static_assert(std::is_const<ParameterType>::value == false,
                "You tried to add a parameter using a type "
                "that is const. Use a non-const type.");

  declare_entry(entry,
                Patterns::Tools::Convert<ParameterType>::to_string(
                  parameter, pattern.clone()),
                pattern,
                documentation,
                has_to_be_set);

  std::string        path = get_current_full_path(entry);
  const unsigned int pattern_index =
    entries->get<unsigned int>(path + path_separator + "pattern");

  auto action = [&, pattern_index](const std::string &val) {
    parameter = Patterns::Tools::Convert<ParameterType>::to_value(
      val, patterns[pattern_index]->clone());
  };
  add_action(entry, action);
}

DEAL_II_NAMESPACE_CLOSE

#endif
