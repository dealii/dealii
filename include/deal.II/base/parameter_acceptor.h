//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#ifndef dealii_base_parameter_acceptor_h
#define dealii_base_parameter_acceptor_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/smartpointer.h>

#include <boost/signals2/signal.hpp>

#include <typeinfo>

DEAL_II_NAMESPACE_OPEN

/**
 * A parameter acceptor base class. This class is used to define a public
 * interface for classes which want to use a single global ParameterHandler to
 * handle parameters. This class declares one static ParameterHandler, and two
 * static functions (declare_all_parameters() and parse_all_parameters()) that
 * manage all of the derived classes.
 *
 * The basic interface provides two subscription mechanisms: a **global
 * subscription mechanism** and a **local subscription mechanism**.
 *
 * The global subscription mechanism is such that whenever an object of a class
 * derived from ParameterAcceptor is created, then a pointer to that
 * object-of-derived-type is registered, together with a path in the parameter
 * file.
 *
 * Such registry is traversed upon invocation of the single function
 * ParameterAcceptor::initialize("file.prm") which in turn calls the method
 * ParameterAcceptor::declare_parameters() for each of the registered classes,
 * reads the file `file.prm` and subsequently calls the method
 * ParameterAcceptor::parse_parameters(), again for each of the registered
 * classes. The method log_info() can be used to extract information about the
 * classes that have been derived from ParameterAcceptor, and that will be
 * parsed when calling ParameterAcceptor::initialize().
 *
 * ParameterAcceptor can be used in three different ways: by overloading the
 * ParameterAcceptor::declare_parameters() and
 * ParameterAcceptor::parse_parameters() methods, by calling its
 * ParameterAcceptor::add_parameter() method for each parameter we want to
 * have, or by constructing a ParameterAcceptorProxy class with your own class,
 * provided that your class implements the @p declare_parameters and
 * @p parse_parameters functions (the first can be a static member in this
 * case).
 *
 * By using the add_parameter method, ParameterAcceptor makes sure that the
 * given parameter is registered in the global parameter handler (by calling
 * ParameterHandler::add_parameter()), at the correct path. If you define all
 * your parameters using the ParameterAcceptor::add_parameter() method, then
 * you don't need to overload any of the virtual methods of this class.
 *
 * If some post processing is required on the parsed values, the user can
 * attach a signal to ParameterAcceptor::declare_parameters_call_back and
 * ParameterAcceptor::parse_parameters_call_back, that are called just after
 * the declare_parameters() and parse_parameters() functions of each derived
 * class. step-69 has an example of doing this.
 *
 * A typical usage of this class is the following:
 *
 * @code
 * // This is your own class, derived from ParameterAcceptor
 * class MyClass : public ParameterAcceptor
 * {
 *   // The constructor of ParameterAcceptor requires a std::string,
 *   // which defines the section name where the parameters of MyClass
 *   // will be stored.
 *   MyClass()
 *     : ParameterAcceptor("Some class name")
 *   {
 *     add_parameter("A param", member_var);
 *   }
 *
 * private:
 *   std::vector<unsigned int> member_var;
 *   ...
 * };
 *
 * int main()
 * {
 *   // Make sure you create your object BEFORE calling
 *   // ParameterAcceptor::initialize()
 *   MyClass class;
 *
 *   // With this call, all derived classes will have their
 *   // parameters initialized
 *   ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 *
 * An implementation that uses user defined declare and parse functions is given
 * by the following example:
 *
 * @code
 * // Again your own class, derived from ParameterAcceptor
 * //
 * // If you don't pass anything to the constructor of
 * // ParameterAcceptor, then the class name is used, "MyClass"
 * // in this case
 * class MyClass : public ParameterAcceptor
 * {
 *   virtual void declare_parameters(ParameterHandler &prm)
 *   {
 *     ...
 *   }
 *
 *   virtual void parse_parameters(ParameterHandler &prm)
 *   {
 *     ...
 *   }
 * };
 *
 * int main()
 * {
 *   // Make sure you create your object BEFORE calling
 *   // ParameterAcceptor::initialize()
 *   MyClass class;
 *   ParameterAcceptor::initialize("file.prm");
 *   class.run();
 * }
 * @endcode
 *
 *
 * Parameter files can be organised into section/subsection/subsubsection.
 * To do so, the std::string passed to ParameterAcceptor within the
 * constructor of the derived class needs to contain the separator "/".
 * In fact, "first/second/third/My Class" will organize the parameters
 * as follows
 *
 * @code
 * subsection first
 *   subsection second
 *     subsection third
 *       subsection My Class
 *        ... # all the parameters
 *       end
 *     end
 *   end
 * end
 * @endcode
 *
 * In the following examples, we propose some use cases with increasing
 * complexities.
 *
 * MyClass is derived from ParameterAcceptor and has a
 * member object that is derived itself from ParameterAcceptor.
 * @code
 * class MyClass : public ParameterAcceptor
 * {
 *   MyClass (std::string name);
 *   virtual void declare_parameters(ParameterHandler &prm);
 * private:
 *   SomeParsedClass<dim> my_subclass;
 *   ...
 * };
 *
 * MyClass::MyClass(std::string name)
 *   : ParameterAcceptor(name)
 *   , my_subclass("Forcing term")
 * {}
 *
 * void MyClass::declare_parmeters(ParameterHandler &prm)
 * {
 *   // many add_parameter(...);
 * }
 *
 * ...
 *
 * int main()
 * {
 *   MyClass mc("My Class");
 *   ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 *
 * In this case, the structure of the parameters will be
 * @code
 * subsection Forcing term
 * ... #parameters of SomeParsedClass
 * end
 * subsection My class
 * ... #all the parameters of MyClass defined in declare_parameters
 * end
 * @endcode
 *
 * Now suppose that in the main file we need two or more objects of MyClass
 * @code
 * int main()
 * {
 *   MyClass ca("Class A");
 *   MyClass cb("Class B");
 *   ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 *
 * What we will read in the parameter file looks like
 * @code
 * subsection Class A
 * ...
 * end
 * subsection Class B
 * ...
 * end
 * subsection Forcing term
 * ...
 * end
 * @endcode
 * Note that there is only one section "Forcing term", this is because
 * both objects have defined the same name for the section of their
 * SomeParsedClass. There are two strategies to change this behaviour. The
 * first one (not recommended) would be to change the name of the section
 * of SomeParsedClass such that it contains also the string passed to
 * the constructor of MyClass:
 * @code
 * MyClass::MyClass(std::string name)
 *  : ParameterAcceptor(name)
 *  , my_subclass(name+" --- forcing term")
 * {}
 * @endcode
 *
 * The other way to proceed (recommended) is to use exploit the
 * /section/subsection approach **in the main class**.
 * @code
 * int main()
 * {
 *   MyClass ca("/Class A/Class");
 *   MyClass cb("/Class B/Class");
 *   ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 * Now, in the parameter file we can find
 * @code
 * subsection Class A
 *   subsection Class
 *   ...
 *   end
 *   subsection Forcing term
 *   ...
 *   end
 * end
 * subsection Class B
 *   subsection Class
 *   ...
 *   end
 *   subsection Forcing term
 *   ...
 *   end
 * end
 * @endcode
 *
 * Note the "/" at the begin of the string name. This is interpreted by
 * ParameterAcceptor like the root folder in Unix systems. The sections "Class
 * A" and "Class B" will not be nested under any section. On the other hand, if
 * the string does not begin with a "/" as in the previous cases the section
 * will be created **under the current path**, which depends on the previously
 * defined sections/subsections/subsubsections. Indeed, the section "Forcing
 * term" is nested under "Class A" or "Class B". To make things more clear.
 * let's consider the following two examples
 *
 * @code
 * int main()
 * {
 *   MyClass ca("/Class A/Class");
 *   MyClass cb("Class B/Class");
 *   ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 * The parameter file will have the following structure
 * @code
 * subsection Class A
 *   subsection Class
 *   ...
 *   end
 *   subsection Forcing term
 *   ...
 *   end
 *   subsection Class B
 *     subsection Class
 *     ...
 *     end
 *     subsection Forcing term
 *     ...
 *     end
 *   end
 * end
 * @endcode
 *
 * If instead one of the paths ends with "/" instead of just
 * a name of the class, subsequent classes will interpret this as
 * a full path, interpreting the class name as a directory name:
 * @code
 * int main()
 * {
 *   MyClass ca("/Class A/Class/");
 *   MyClass cb("Class B/Class");
 *   ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 * The parameter file will have the following structure
 * @code
 * subsection Class A
 *   subsection Class
 *      ...
 *      subsection Forcing term
 *      ...
 *      end
 *      subsection Class B
 *          subsection Class
 *          ...
 *          end
 *          subsection Forcing term
 *          ...
 *          end
 *      end
 *   end
 * end
 * @endcode
 *
 * As a final remark, in order to allow a proper management of all the
 * sections/subsections, the instantiation of objects and the call to
 * ParameterAcceptor::initialize() cannot be done on multiple, concurrently
 * running threads.
 *
 * If you pass an empty name, the boost::core::demangle() function is used to
 * fill the section name with a human readable version of the class name
 * itself.
 *
 * See the tutorial program step-60 for an example on how to use this class.
 *
 * @author Luca Heltai, 2017.
 */
class ParameterAcceptor : public Subscriptor
{
public:
  /**
   * The constructor adds derived classes to the list of acceptors. If
   * a section name is specified, then this is used to scope the
   * parameters in the given section, otherwise a pretty printed
   * version of the derived class is used.
   */
  ParameterAcceptor(const std::string &section_name = "");

  /**
   * Destructor.
   */
  virtual ~ParameterAcceptor() override;

  /**
   * Call declare_all_parameters(), read the parameters from `filename` (only
   * if `filename` is a non-empty string), and then call
   * parse_all_parameters().
   *
   * If the parameter `filename` is the empty string, then no attempt to read a
   * parameter file is done. This may be useful if you are ok with using
   * default values, and don't want to read external files to use a class
   * derived from ParameterAcceptor.
   *
   * If @p output_filename is not the empty string, then we write the content
   * that was read into the @p output_filename file, using the style specified
   * in @p output_style_for_output_filename. The format of both input and output
   * files are selected using the extensions of the files themselves. This can
   * be either `prm`, `xml`, or `json` for the @p filename, and any of the
   * supported formats for the @p output_filename.
   *
   * If the input file does not exist, a default one with the same name is
   * created for you following the style specified in
   * @p output_style_for_filename, and an exception is thrown.
   *
   * By default, the file format used to write the files is deduced from
   * the extension of the file names. If the corresponding
   * ParameterHandler::OutputStyle specifies a format specification, this must
   * be compatible with the file extension, or an exception will be thrown.
   *
   * If the extension is not recognized, and you do not specify a format in the
   * corresponding ParameterHandler::OutputStyle, an assertion is thrown.
   *
   * @param filename Input file name
   * @param output_filename Output file name
   * @param output_style_for_output_filename How to write the output file
   * @param prm The ParameterHandler to use
   * @param output_style_for_filename How to write the default input file if it
   * does not exist
   */
  static void
  initialize(const std::string &filename        = "",
             const std::string &output_filename = "",
             const ParameterHandler::OutputStyle
                                                 output_style_for_output_filename = ParameterHandler::Short,
             ParameterHandler &                  prm = ParameterAcceptor::prm,
             const ParameterHandler::OutputStyle output_style_for_filename =
               ParameterHandler::DefaultStyle);

  /**
   * Call declare_all_parameters(), read the parameters from the `input_stream`
   * in `prm` format, and then call parse_all_parameters().
   *
   * An exception is thrown if the `input_stream` is invalid.
   *
   * @param input_stream Input stream
   * @param prm The ParameterHandler to use
   */
  static void
  initialize(std::istream &    input_stream,
             ParameterHandler &prm = ParameterAcceptor::prm);


  /**
   * Clear class list and global parameter file.
   */
  static void
  clear();

  /**
   * Derived classes can use this method to declare their parameters.
   * ParameterAcceptor::initialize() calls it for each derived class. The
   * default implementation is empty.
   */
  virtual void
  declare_parameters(ParameterHandler &prm);

  /**
   * Declare parameter call back. This signal is triggered right after
   * declare_parameters() has been called, to allow users to prepare their
   * variables right after parameters have been decalred. The default
   * implementation is empty.
   */
  boost::signals2::signal<void()> declare_parameters_call_back;

  /**
   * Derived classes can use this method to parse their parameters.
   * ParameterAcceptor::initialize() calls it for each derived class. The
   * default implementation is empty.
   */
  virtual void
  parse_parameters(ParameterHandler &prm);

  /**
   * Parse parameter call back. This function is called at the end of
   * parse_parameters(), to allow users to process their parameters right after
   * they have been parsed. The default implementation is empty.
   *
   * You can use this function, for example, to create a quadrature rule after
   * you have read how many quadrature points you wanted to use from the
   * parameter file.
   */
  boost::signals2::signal<void()> parse_parameters_call_back;

  /**
   * Parse the given ParameterHandler. This function enters the
   * subsection returned by get_section_name() for each derived class,
   * and parses all parameters that were added using add_parameter().
   */
  static void
  parse_all_parameters(ParameterHandler &prm = ParameterAcceptor::prm);

  /**
   * Initialize the global ParameterHandler with all derived classes
   * parameters.This function enters the subsection returned by
   * get_section_name() for each derived class, and declares all parameters
   * that were added using add_parameter().
   */
  static void
  declare_all_parameters(ParameterHandler &prm = ParameterAcceptor::prm);

  /**
   * Return the section name of this class. If a name was provided
   * at construction time, then that name is returned, otherwise it
   * returns the demangled name of this class.
   */
  std::string
  get_section_name() const;

  /**
   * Traverse all registered classes, and figure out what subsections we need to
   * enter.
   */
  std::vector<std::string>
  get_section_path() const;

  /**
   * Add a parameter in the correct path. This method forwards all arguments to
   * the prm.add_parameter() method, after entering the correct section path.
   * By default it uses the ParameterAcceptor::prm variable as
   * ParameterHandler.
   *
   * See the documentation of ParameterHandler::add_parameter() for more
   * information.
   */
  template <class ParameterType>
  void
  add_parameter(const std::string &          entry,
                ParameterType &              parameter,
                const std::string &          documentation = "",
                ParameterHandler &           prm_          = prm,
                const Patterns::PatternBase &pattern =
                  *Patterns::Tools::Convert<ParameterType>::to_pattern());

  /**
   * The global parameter handler.
   */
  static ParameterHandler prm;

  /**
   * Make sure we enter the right subsection of the given parameter.
   */
  void
  enter_my_subsection(ParameterHandler &prm);

  /**
   * This function undoes what the enter_my_subsection() function did. It only
   * makes sense if enter_my_subsection() was called on `prm` before this one.
   */
  void
  leave_my_subsection(ParameterHandler &prm);

private:
  /**
   * A list containing all constructed classes of type
   * ParameterAcceptor.
   */
  static std::vector<SmartPointer<ParameterAcceptor>> class_list;

  /** The index of this specific class within the class list. */
  const unsigned int acceptor_id;

  /**
   * Separator between sections.
   */
  static const char sep = '/';

protected:
  /** The subsection name for this class. */
  const std::string section_name;
};



/**
 * A proxy ParameterAcceptor wrapper for classes that have a static member
 * function @p declare_parameters, and a non virtual @p parse_parameters method.
 *
 * If you cannot or do not want to derive your "parameter accepting" class from
 * ParameterAcceptor, for example if by design you are required to have a
 * static member function @p declare_parameters and a member @p
 * parse_parameters, or if someone has already implemented such a class for
 * you, and only provides you with an API that you cannot modify, then you may
 * be able to use ParameterAcceptor facilities nonetheless, by wrapping your
 * class into ParameterAcceptorProxy.
 *
 * This class implements the public interface of ParameterAcceptor, and at the
 * same time it derives from the template class @p SourceClass, allowing you to
 * register your existing @p SourceClass as a ParameterAcceptor class, without
 * requiring you to explicitly derive your @p SourceClass from
 * ParameterAcceptor.
 *
 * An example usage is given by the following snippet of code, using
 * Functions::ParsedFunction as an example source class:
 *
 * @code
 * ParameterAcceptorProxy<Functions::ParsedFunction<2> > fun("Some function");
 * ParameterAcceptor::initialize("test.prm");
 * @endcode
 *
 * The above snippet of code will initialize ParameterAcceptor::prm with a
 * section "Some function", and will correctly parse and assign to the object
 * `fun` the expression parsed from the file `test.prm`. If non-existent, the
 * program will exit, and generate it for you (here you can see the resulting
 * short text version of the parameter file generated with the above snippet):
 *
 * @code
 * # Parameter file generated with
 * # DEAL_II_PACKAGE_VERSION = 9.0.0-pre
 * subsection Some function
 *   set Function constants  =
 *   set Function expression = 0
 *   set Variable names      = x,y,t
 * end
 * @endcode
 *
 * The resulting `fun` object, is both a ParsedFunction object and a
 * ParameterAcceptor one, allowing you to use it as a replacement of the
 * ParsedFunction class, with automatic declaration and parsing of parameter
 * files.
 *
 * See the tutorial program step-60 for an example on how to use this class.
 *
 * @author Luca Heltai, 2018
 */
template <class SourceClass>
class ParameterAcceptorProxy : public SourceClass, public ParameterAcceptor
{
public:
  /**
   * Default constructor. The argument `section_name` is forwarded to the
   * constructor of the ParameterAcceptor class, while all other arguments
   * are passed to the SourceClass constructor.
   */
  template <typename... Args>
  ParameterAcceptorProxy(const std::string &section_name, Args... args);

  /**
   * Overloads the ParameterAcceptor::declare_parameters function, by calling
   * @p SourceClass::declare_parameters with @p prm as an argument.
   */
  virtual void
  declare_parameters(ParameterHandler &prm) override;

  /**
   * Overloads the ParameterAcceptor::parse_parameters function, by calling
   * @p SourceClass::parse_parameters with @p prm as an argument.
   */
  virtual void
  parse_parameters(ParameterHandler &prm) override;
};



// Inline and template functions
template <class ParameterType>
void
ParameterAcceptor::add_parameter(const std::string &          entry,
                                 ParameterType &              parameter,
                                 const std::string &          documentation,
                                 ParameterHandler &           prm,
                                 const Patterns::PatternBase &pattern)
{
  enter_my_subsection(prm);
  prm.add_parameter(entry, parameter, documentation, pattern);
  leave_my_subsection(prm);
}



template <class SourceClass>
template <typename... Args>
ParameterAcceptorProxy<SourceClass>::ParameterAcceptorProxy(
  const std::string &section_name,
  Args... args)
  : SourceClass(args...)
  , ParameterAcceptor(section_name)
{}



template <class SourceClass>
void
ParameterAcceptorProxy<SourceClass>::declare_parameters(ParameterHandler &prm)
{
  SourceClass::declare_parameters(prm);
}



template <class SourceClass>
void
ParameterAcceptorProxy<SourceClass>::parse_parameters(ParameterHandler &prm)
{
  SourceClass::parse_parameters(prm);
}

DEAL_II_NAMESPACE_CLOSE

#endif
