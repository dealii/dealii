// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
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

#ifndef dealii_patterns_h
#define dealii_patterns_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx14/algorithm.h>
#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/base/std_cxx14/utility.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/component_mask.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/archive/basic_archive.hpp>
#include <boost/core/demangle.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/property_tree/ptree_serialization.hpp>
#include <boost/serialization/split_member.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <array>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations for interfaces and friendship
#ifndef DOXYGEN
class LogStream;
class MultipleParameterLoop;
template <int dim>
class FunctionParser;
#endif

/**
 * Namespace for a few classes that act as patterns for the ParameterHandler
 * class. These classes implement an interface that checks whether a parameter
 * in an input file matches a certain pattern, such as "being boolean", "an
 * integer value", etc.
 *
 * @ingroup input
 */
namespace Patterns
{
  /**
   * Base class to declare common interface. The purpose of this class is
   * mostly to define the interface of patterns, and to force derived classes
   * to have a <tt>clone</tt> function. It is thus, in the languages of the
   * "Design Patterns" book (Gamma et al.), a "prototype".
   */
  class PatternBase
  {
  public:
    /**
     * Make destructor of this and all derived classes virtual.
     */
    virtual ~PatternBase() = default;

    /**
     * Return <tt>true</tt> if the given string matches the pattern.
     */
    virtual bool
    match(const std::string &test_string) const = 0;

    /**
     * List of possible description output formats.
     *
     * Capitalization chosen for similarity to ParameterHandler::OutputStyle.
     */
    enum OutputStyle
    {
      /**
       * Simple text suitable for machine parsing in the static public member
       * functions for all of the built in inheriting classes.
       *
       * Preferably human readable, but machine parsing is more critical.
       */
      Machine,
      /**
       * Easily human readable plain text format suitable for plain text
       * documentation.
       */
      Text,
      /**
       * Easily human readable LaTeX format suitable for printing in manuals.
       */
      LaTeX
    };

    /**
     * Return a string describing the pattern.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const = 0;

    /**
     * Return a pointer to an exact copy of the object. This is necessary
     * since we want to store objects of this type in containers, were we need
     * to copy objects without knowledge of their actual data type (we only
     * have pointers to the base class).
     *
     * Ownership of the objects returned by this function is passed to the
     * caller of this function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const = 0;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object. To avoid unnecessary overhead, we do not force derived classes
     * to provide this function as a virtual overloaded one, but rather try to
     * cast the present object to one of the known derived classes and if that
     * fails then take the size of this base class instead and add 32 byte
     * (this value is arbitrary, it should account for virtual function
     * tables, and some possible data elements). Since there are usually not
     * many thousands of objects of this type around, and since the
     * memory_consumption mechanism is used to find out where memory in the
     * range of many megabytes is, this seems like a reasonable approximation.
     *
     * On the other hand, if you know that your class deviates from this
     * assumption significantly, you can still overload this function.
     */
    virtual std::size_t
    memory_consumption() const;
  };

  /**
   * Return pointer to the correct derived class based on description.
   */
  std::unique_ptr<PatternBase>
  pattern_factory(const std::string &description);

  namespace internal
  {
    /**
     * Escape the string @p input for the specified @p style so that characters
     * will appear as intended. For example, characters like _ can not be
     * written as is in LateX and have to be escaped as \_.
     */
    std::string
    escape(const std::string &input, const PatternBase::OutputStyle style);

  } // namespace internal

  /**
   * Test for the string being an integer. If bounds are given to the
   * constructor, then the integer given also needs to be within the interval
   * specified by these bounds. Note that unlike common convention in the C++
   * standard library, both bounds of this interval are inclusive; the reason
   * is that in practice in most cases, one needs closed intervals, but these
   * can only be realized with inclusive bounds for non-integer values. We
   * thus stay consistent by always using closed intervals.
   *
   * If the upper bound given to the constructor is smaller than the
   * lower bound, then every integer is allowed.
   *
   * Giving bounds may be useful if for example a value can only be positive
   * and less than a reasonable upper bound (for example the number of
   * refinement steps to be performed), or in many other cases.
   */
  class Integer : public PatternBase
  {
  public:
    /**
     * Minimal integer value. If the numeric_limits class is available use
     * this information to obtain the extremal values, otherwise set it so
     * that this class understands that all values are allowed.
     */
    static const int min_int_value;

    /**
     * Maximal integer value. If the numeric_limits class is available use
     * this information to obtain the extremal values, otherwise set it so
     * that this class understands that all values are allowed.
     */
    static const int max_int_value;

    /**
     * Constructor. Bounds can be specified within which a valid
     * parameter has to be. If the upper bound is smaller than the
     * lower bound, then the entire set of integers is implied. The
     * default values are chosen such that no bounds are enforced on
     * parameters.
     *
     * Note that the range implied by an object of the current type
     * is inclusive of both bounds values, i.e., the @p upper_bound is
     * an allowed value, rather than indicating a half-open value as
     * is often done in other contexts.
     */
    Integer(const int lower_bound = min_int_value,
            const int upper_bound = max_int_value);

    /**
     * Return <tt>true</tt> if the string is an integer and its value is
     * within the specified range.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. If bounds were specified to the constructor, then include them
     * into this description.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Integer>
    create(const std::string &description);

  private:
    /**
     * Value of the lower bound. A number that satisfies the
     * @ref match
     * operation of this class must be equal to this value or larger, if the
     * bounds of the interval for a valid range.
     */
    const int lower_bound;

    /**
     * Value of the upper bound. A number that satisfies the
     * @ref match
     * operation of this class must be equal to this value or less, if the
     * bounds of the interval for a valid range.
     */
    const int upper_bound;

    /**
     * Initial part of description
     */
    static const char *description_init;
  };

  /**
   * Test for the string being a <tt>double</tt>. If bounds are given to the
   * constructor, then the integer given also needs to be within the interval
   * specified by these bounds. Note that unlike common convention in the C++
   * standard library, both bounds of this interval are inclusive; the reason
   * is that in practice in most cases, one needs closed intervals, but these
   * can only be realized with inclusive bounds for non-integer values. We
   * thus stay consistent by always using closed intervals.
   *
   * If the upper bound given to the constructor is smaller than the
   * lower bound, then every double precision number is allowed.
   *
   * Giving bounds may be useful if for example a value can only be positive
   * and less than a reasonable upper bound (for example damping parameters
   * are frequently only reasonable if between zero and one), or in many other
   * cases.
   */
  class Double : public PatternBase
  {
  public:
    /**
     * Minimal double value used as default value, taken from
     * <tt>std::numeric_limits</tt>.
     */
    static const double min_double_value;

    /**
     * Maximal double value used as default value, taken from
     * <tt>std::numeric_limits</tt>.
     */
    static const double max_double_value;

    /**
     * Constructor. Bounds can be specified within which a valid
     * parameter has to be. If the upper bound is smaller than the
     * lower bound, then the entire set of double precision numbers is
     * implied. The default values are chosen such that no bounds are
     * enforced on parameters.
     */
    Double(const double lower_bound = min_double_value,
           const double upper_bound = max_double_value);

    /**
     * Return <tt>true</tt> if the string is a number and its value is within
     * the specified range.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. If bounds were specified to the constructor, then include them
     * into this description.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Creates a new object on the heap using @p new if the given
     * @p description is a valid format (for example created by calling
     * description() on an existing object), or @p nullptr otherwise. Ownership
     * of the returned object is transferred to the caller of this function,
     * which should be freed using @p delete.
     */
    static std::unique_ptr<Double>
    create(const std::string &description);

  private:
    /**
     * Value of the lower bound. A number that satisfies the
     * @ref match
     * operation of this class must be equal to this value or larger, if the
     * bounds of the interval form a valid range.
     */
    const double lower_bound;

    /**
     * Value of the upper bound. A number that satisfies the
     * @ref match
     * operation of this class must be equal to this value or less, if the
     * bounds of the interval form a valid range.
     */
    const double upper_bound;

    /**
     * Initial part of description
     */
    static const char *description_init;
  };

  /**
   * Test for the string being one of a sequence of values given like a
   * regular expression. For example, if the string given to the constructor
   * is <tt>"red|blue|black"</tt>, then the
   * @ref match
   * function returns <tt>true</tt> exactly if the string is either "red" or
   * "blue" or "black". Spaces around the pipe signs do not matter and are
   * eliminated.
   */
  class Selection : public PatternBase
  {
  public:
    /**
     * Constructor. Take the given parameter as the specification of valid
     * strings.
     */
    Selection(const std::string &seq);

    /**
     * Return <tt>true</tt> if the string is an element of the description
     * list passed to the constructor.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the list of valid strings passed to the
     * constructor.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const override;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Selection>
    create(const std::string &description);

  private:
    /**
     * List of valid strings as passed to the constructor. We don't make this
     * string constant, as we process it somewhat in the constructor.
     */
    std::string sequence;

    /**
     * Initial part of description
     */
    static const char *description_init;
  };


  /**
   * This pattern matches a list of values separated by commas (or another
   * string), each of which have to match a pattern given to the constructor.
   * With two additional parameters, the number of elements this list has to
   * have can be specified. If none is specified, the list may have zero or
   * more entries.
   */
  class List : public PatternBase
  {
  public:
    /**
     * Maximal integer value. If the numeric_limits class is available use
     * this information to obtain the extremal values, otherwise set it so
     * that this class understands that all values are allowed.
     */
    static const unsigned int max_int_value;

    /**
     * Constructor. Take the given parameter as the specification of valid
     * elements of the list.
     *
     * The three other arguments can be used to denote minimal and maximal
     * allowable lengths of the list, and the string that is used as a
     * separator between elements of the list.
     */
    List(const PatternBase &base_pattern,
         const unsigned int min_elements = 0,
         const unsigned int max_elements = max_int_value,
         const std::string &separator    = ",");


    /**
     * Return the internally stored separator.
     */
    const std::string &
    get_separator() const;

    /**
     * Return the internally stored base pattern.
     */
    const PatternBase &
    get_base_pattern() const;

    /**
     * Copy constructor.
     */
    List(const List &other);

    /**
     * Return <tt>true</tt> if the string is a comma-separated list of strings
     * each of which match the pattern given to the constructor.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<List>
    create(const std::string &description);

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const override;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception.
     */
    DeclException2(ExcInvalidRange,
                   int,
                   int,
                   << "The values " << arg1 << " and " << arg2
                   << " do not form a valid range.");
    //@}
  private:
    /**
     * Copy of the pattern that each element of the list has to satisfy.
     */
    std::unique_ptr<PatternBase> pattern;

    /**
     * Minimum number of elements the list must have.
     */
    const unsigned int min_elements;

    /**
     * Maximum number of elements the list must have.
     */
    const unsigned int max_elements;

    /**
     * Separator between elements of the list.
     */
    const std::string separator;

    /**
     * Initial part of description
     */
    static const char *description_init;
  };


  /**
   * This pattern matches a list of comma-separated values each of which
   * denotes a pair of key and value. Both key and value have to match a
   * pattern given to the constructor. For each entry of the map, parameters
   * have to be entered in the form <code>key: value</code>. In other words, a
   * map is described in the form <code>key1: value1, key2: value2, key3:
   * value3, ...</code>. Two constructor arguments allow to choose a delimiter
   * between pairs other than the comma, and a delimiter between key and value
   * other than colon.
   *
   * With two additional parameters, the number of elements this list has to
   * have can be specified. If none is specified, the map may have zero or
   * more entries.
   */
  class Map : public PatternBase
  {
  public:
    /**
     * Maximal integer value. If the numeric_limits class is available use
     * this information to obtain the extremal values, otherwise set it so
     * that this class understands that all values are allowed.
     */
    static const unsigned int max_int_value;

    /**
     * Constructor. Take the given parameter as the specification of valid
     * elements of the list.
     *
     * The four other arguments can be used to denote minimal and maximal
     * allowable lengths of the list as well as the separators used to delimit
     * pairs of the map and the symbol used to separate keys and values.
     */
    Map(const PatternBase &key_pattern,
        const PatternBase &value_pattern,
        const unsigned int min_elements        = 0,
        const unsigned int max_elements        = max_int_value,
        const std::string &separator           = ",",
        const std::string &key_value_separator = ":");

    /**
     * Copy constructor.
     */
    Map(const Map &other);

    /**
     * Return <tt>true</tt> if the string is a comma-separated list of strings
     * each of which match the pattern given to the constructor.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Map>
    create(const std::string &description);

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const override;

    /**
     * Return a reference to the key pattern.
     */
    const PatternBase &
    get_key_pattern() const;

    /**
     * Return a reference to the value pattern.
     */
    const PatternBase &
    get_value_pattern() const;

    /**
     * Return the separator of the map entries.
     */
    const std::string &
    get_separator() const;

    /**
     * Return the key-value separator.
     */
    const std::string &
    get_key_value_separator() const;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception.
     */
    DeclException2(ExcInvalidRange,
                   int,
                   int,
                   << "The values " << arg1 << " and " << arg2
                   << " do not form a valid range.");
    //@}
  private:
    /**
     * Copy of the patterns that each key and each value of the map has to
     * satisfy.
     */
    std::unique_ptr<PatternBase> key_pattern;
    std::unique_ptr<PatternBase> value_pattern;

    /**
     * Minimum number of elements the list must have.
     */
    const unsigned int min_elements;

    /**
     * Maximum number of elements the list must have.
     */
    const unsigned int max_elements;

    /**
     * Separator between elements of the list.
     */
    const std::string separator;


    /**
     * Separator between keys and values.
     */
    const std::string key_value_separator;

    /**
     * Initial part of description
     */
    static const char *description_init;
  };



  /**
   * This pattern matches colon-separated values of arbitrary types. Each type
   * has to match a pattern given to the constructor.
   *
   * An example usage is the following:
   *
   * @code
   * std::vector< std::unique_ptr<Patterns::PatternBase> > ps;
   *
   * ps.push_back(std::unique_ptr<Patterns::Integer>());
   * ps.push_back(std::unique_ptr<Patterns::Double>());
   * ps.push_back(std::unique_ptr<Patterns::Anything>());
   *
   * Patterns::Tuple pattern(ps, ":");
   *
   * bool check = ps.match("5 : 3.14 : Ciao"); // check = true
   * @endcode
   *
   * or, if you want to exploit ParameterHandler::add_parameter():
   *
   * @code
   * using T = std::tuple<std::string, Point<3>, unsigned int>;
   *
   * T a = Patterns::Tools::Convert<T>::to_value("Ciao : 1.0, 2.0, 3.0 : 33");
   *
   * ParameterHandler prm;
   * prm.add_parameter("A tuple", a);
   *
   * prm.log_parameters(deallog);
   * // DEAL:parameters::A tuple: Ciao : 1.000000, 2.000000, 3.000000 : 33
   *
   * prm.set("A tuple", "Mondo : 2.0, 3.0, 4.0 : 34");
   * prm.log_parameters(deallog);
   * // DEAL:parameters::A tuple: Mondo : 2.0, 3.0, 4.0 : 34
   *
   * deallog << Patterns::Tools::Convert<T>::to_string(a) << std::endl;
   * // DEAL::Mondo : 2.000000, 3.000000, 4.000000 : 34
   * @endcode
   *
   * The constructor expects a vector of Patterns, and optionally a string
   * specifying the separator to use when parsing the Tuple from a string.
   *
   * The default separator is a colon, owing to the fact that a pair is in fact
   * a tuple with two elements.
   *
   * @author Luca Heltai, 2017.
   */
  class Tuple : public PatternBase
  {
  public:
    /**
     * Constructor. Use a vector of unique pointers to Patterns to construct
     * the tuple.
     *
     * @param patterns The pattern each object of the Tuple should match
     * @param separator An optional string used to delimit each element
     * Constructor.
     */
    Tuple(const std::vector<std::unique_ptr<PatternBase>> &patterns,
          const std::string &                              separator = ":");

    /**
     * Constructor. Same as above, specialized for const char *. This is
     * necessary to avoid compilers errors due to the variadic constructors
     * provided below.
     */
    Tuple(const std::vector<std::unique_ptr<PatternBase>> &patterns,
          const char *                                     separator);


    /**
     * Constructor. Creates a Tuple from more than one class derived from
     * PatternBase.
     *
     * @param separator What separator to use.
     * @param patterns The list of patterns to use
     */
    template <class... PatternTypes>
    Tuple(const std::string &separator, const PatternTypes &... patterns);

    /**
     * Constructor. This is needed to allow users to specify
     * directly the separator without using std::string(";").
     *
     * Since we support a pure variadic templates version, without this
     * specialization, the compiler will fail with cryptic errors.
     */
    template <class... PatternTypes>
    Tuple(const char *separator, const PatternTypes &... patterns);

    /**
     * Constructor. Same as above, using the default separator.
     *
     * @param patterns The list of patterns to use
     */
    template <typename... Patterns>
    Tuple(const Patterns &... patterns);

    /**
     * Copy constructor.
     */
    Tuple(const Tuple &other);

    /**
     * Return <tt>true</tt> if the string is a list of strings
     * each of which matches the patterns given to the constructor.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Tuple>
    create(const std::string &description);

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const override;

    /**
     * Return a reference to the i-th pattern in the tuple.
     */
    const PatternBase &
    get_pattern(const unsigned int i) const;

    /**
     * Return the separator of the tuple entries.
     */
    const std::string &
    get_separator() const;

  private:
    /**
     * Copy of the patterns stored in the Tuple.
     */
    std::vector<std::unique_ptr<PatternBase>> patterns;

    /**
     * Separator between elements of the list.
     */
    const std::string separator;

    /**
     * Initial part of description.
     */
    static const char *description_init;
  };


  /**
   * This class is much like the Selection class, but it allows the input to
   * be a comma-separated list of values which each have to be given in the
   * constructor argument. The input is allowed to be empty or contain values
   * more than once and have an arbitrary number of spaces around commas. Of
   * course commas are not allowed inside the values given to the constructor.
   *
   * For example, if the string to the constructor was <tt>"ucd|gmv|eps"</tt>,
   * then the following would be legal inputs: "eps", "gmv, eps", or "".
   */
  class MultipleSelection : public PatternBase
  {
  public:
    /**
     * Constructor. @p seq is a list of valid options separated by "|".
     */
    MultipleSelection(const std::string &seq);

    /**
     * Return <tt>true</tt> if the string is an element of the description
     * list passed to the constructor.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the list of valid strings passed to the
     * constructor.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<MultipleSelection>
    create(const std::string &description);

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const override;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception.
     */
    DeclException1(
      ExcCommasNotAllowed,
      int,
      << "A comma was found at position " << arg1
      << " of your input string, but commas are not allowed here.");
    //@}
  private:
    /**
     * List of valid strings as passed to the constructor. We don't make this
     * string constant, as we process it somewhat in the constructor.
     */
    std::string sequence;

    /**
     * Initial part of description
     */
    static const char *description_init;
  };

  /**
   * Test for the string being either "true" or "false". This is mapped to the
   * Selection class.
   */
  class Bool : public Selection
  {
  public:
    /**
     * Constructor.
     */
    Bool();

    /**
     * Return a description of the pattern that valid strings are expected to
     * match.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Bool>
    create(const std::string &description);

  private:
    /**
     * Initial part of description
     */
    static const char *description_init;
  };

  /**
   * Always returns <tt>true</tt> when testing a string.
   */
  class Anything : public PatternBase
  {
  public:
    /**
     * Constructor. (Allow for at least one non-virtual function in this
     * class, as otherwise sometimes no virtual table is emitted.)
     */
    Anything() = default;

    /**
     * Return <tt>true</tt> if the string matches its constraints, i.e.
     * always.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the string <tt>"[Anything]"</tt>.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Anything>
    create(const std::string &description);

  private:
    /**
     * Initial part of description
     */
    static const char *description_init;
  };


  /**
   * A pattern that can be used to indicate when a parameter is intended to be
   * the name of a file. By itself, this class does not check whether the
   * string that is given in a parameter file actually corresponds to an
   * existing file (it could, for example, be the name of a file to which you
   * want to write output). Functionally, the class is therefore equivalent to
   * the Anything class. However, it allows to specify the <i>intent</i> of a
   * parameter. The flag given to the constructor also allows to specify
   * whether the file is supposed to be an input or output file.
   *
   * The reason for the existence of this class is to support graphical user
   * interfaces for editing parameter files. These may open a file selection
   * dialog if the filename is supposed to represent an input file.
   */
  class FileName : public PatternBase
  {
  public:
    /**
     * Files can be used for input or output. This can be specified in the
     * constructor by choosing the flag <tt>type</tt>.
     */
    enum FileType
    {
      /**
       * Open for input.
       */
      input = 0,
      /**
       * Open for output.
       */
      output = 1
    };

    /**
     * Constructor.  The type of the file can be specified by choosing the
     * flag.
     */
    FileName(const FileType type = input);

    /**
     * Return <tt>true</tt> if the string matches its constraints, i.e.
     * always.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the string <tt>"[Filename]"</tt>.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * file type flag
     */
    FileType file_type;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<FileName>
    create(const std::string &description);

  private:
    /**
     * Initial part of description
     */
    static const char *description_init;
  };


  /**
   * A pattern that can be used to indicate when a parameter is intended to be
   * the name of a directory. By itself, this class does not check whether the
   * string that is given in a parameter file actually corresponds to an
   * existing directory. Functionally, the class is therefore equivalent to
   * the Anything class. However, it allows to specify the <i>intent</i> of a
   * parameter.
   *
   * The reason for the existence of this class is to support graphical user
   * interfaces for editing parameter files. These may open a file selection
   * dialog to select or create a directory.
   */
  class DirectoryName : public PatternBase
  {
  public:
    /**
     * Constructor.
     */
    DirectoryName() = default;

    /**
     * Return <tt>true</tt> if the string matches its constraints, i.e.
     * always.
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the string <tt>"[Filename]"</tt>.
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * Create a new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<DirectoryName>
    create(const std::string &description);

  private:
    /**
     * Initial part of description
     */
    static const char *description_init;
  };


  /**
   * Namespace for a few classes and functions that act on values and patterns,
   * and allow to convert from non elementary types to strings and vice versa.
   *
   * A typical usage of these tools is in the following example:
   *
   * @code
   * using T = std::vector<unsigned int>;
   *
   * T vec(3);
   * vec[0] = 1;
   * vec[1] = 3;
   * vec[2] = 5;
   *
   * auto pattern = Patterns::Tools::Convert<T>::to_pattern();
   *
   * std::cout << pattern->description() << std::endl;
   * // [List of <[Integer]> of length 0...4294967295 (inclusive)]
   *
   * auto s = Patterns::Tools::Convert<T>::to_string(vec);
   * std::cout << s << std::endl;
   * // 1, 2, 3
   *
   * auto vec = Patterns::Tools::Convert<T>::to_value("2,3,4,5");
   * // now vec has size 4, and contains the elements 2,3,4,5
   *
   * std::cout << internal::RankInfo<T>::list_rank << std::endl; // Outputs 1
   * std::cout << internal::RankInfo<T>::map_rank  << std::endl; // Outputs 0
   * @endcode
   *
   * Convert<T> is used by the function Patterns::Tools::add_parameter() in this
   * namespace. Internally it uses the internal::RankInfo<T> class to decide how
   * many different separators are required to convert the given type to a
   * string.
   *
   * For example, to write vectors of vectors, the default is to use "," for the
   * first (inner) separator, and ";" for the second (outer) separator, i.e.
   *
   * @code
   * std::vector<std::vector<unsigned int>> vec;
   * vec = Convert<decltype(vec)>::to_value("1,2,3 ; 4,5,6");
   *
   * s = convert<decltype(vec[0])>::to_string(vec[0]);
   * // s now contains the string "1,2,3"
   * @endcode
   *
   * Separators for Patterns::List and Patterns::Map compatible types are
   * selected according to the
   * rank of the list and map objects, using the arrays
   * Patterns::Tools::internal::default_list_separator and
   * Patterns::Tools::internal::default_map_separator.
   *
   * They are currently set to:
   *
   * @code
   * default_list_separator{{","  ,  ";"  ,  "|"  ,   "%"}};
   * default_map_separator {{":"  ,  "="  ,  "@"  ,   "#"}};
   * @endcode
   *
   * When one needs a mixture of Patterns::List and Patterns::Map types, their
   * RankInfo is computed by taking the maximum of the vector_rank of the Key
   * and of the Value type, so that, for example, it is possible to have the
   * following
   * @code
   * ... // Build compare class
   * std::map<std::vector<unsigned int>, std::vector<double>, compare> map;
   *
   * map = convert<decltype(map)>::to_value(
   *   "1,2,3 : 5.0,6.0,7.0  ; 8,9,10 : 11.0,12.0,13.0");
   *
   * @endcode
   *
   * Some non elementary types are supported, like Point(), or
   * std::complex<double>. If you wish to support more types, you have to
   * specialize the Convert struct as well as the RankInfo struct.
   *
   * @ingroup input
   * @author Luca Heltai, 2017
   */
  namespace Tools
  {
    /**
     * Converter class. This class is used to generate strings and Patterns
     * associated to the given type, and to convert from a string to the given
     * type and vice versa.
     *
     * The second template parameter is used internally to allow for advanced
     * SFINAE (substitution failure is not an error) tricks used to specialise
     * this class for arbitrary STL containers and maps.
     *
     * @author Luca Heltai, 2017
     */
    template <class T, class Enable = void>
    struct Convert
    {
      /**
       * Return a std::unique_ptr to a Pattern that can be used to interpret a
       * string as the type of the template argument, and the other way around.
       *
       * While the current function (in the general Convert template) is
       * deleted, it is implemented and available in the specializations of the
       * Convert
       * class template for particular kinds of template arguments @p T.
       */
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern() = delete;

      /**
       * Return a string containing a textual version of the variable s. Use the
       * pattern passed to perform the conversion, or create and use a default
       * one.
       *
       * While the current function (in the general Convert template) is
       * deleted, it is implemented and available in the specializations of the
       * Convert
       * class template for particular kinds of template arguments @p T.
       */
      static std::string
      to_string(const T &                                     s,
                const std::unique_ptr<Patterns::PatternBase> &p =
                  Convert<T>::to_pattern()) = delete;

      /**
       * Convert a string to a value, using the given pattern. Use the pattern
       * passed to perform the conversion, or create and use a default one.
       *
       * While the current function (in the general Convert template) is
       * deleted, it is implemented and available in the specializations of the
       * Convert
       * class template for particular kinds of template arguments @p T.
       */
      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &p =
                 Convert<T>::to_pattern()) = delete;
    };

    /**
     * A utility function that simplifies the conversion to strings of
     * arbitrarily complex types.
     *
     * This function calls the method Convert<T>::to_string() with the default
     * pattern. An example usage is the following:
     *
     * @code
     * auto t = std::make_tuple(1.0, std::make_pair(1, "ciao"));
     * auto s = Patterns::Tools::to_string(t);
     *
     * std::cout << s; // will print "1 % 1 : ciao""
     * @endcode
     *
     * See the documentation of the class Patterns::Tools::Convert, and of the
     * helper class Patterns::Tools::RankInfo for details on the way separators
     * are selected when outputting STL container types.
     *
     * @author Luca Heltai, 2018
     */
    template <typename T>
    std::string
    to_string(const T &t);

    /**
     * A utility function that simplifies the conversion from strings to
     * arbitrary types.
     *
     * This function calls the method Convert<T>::to_value() with the default
     * pattern. An example usage is the following:
     *
     * @code
     * auto t = std::make_tuple(1.0, std::make_pair(1, "ciao"));
     * // replace the value of 't' by the parsed content of the string argument:
     * Patterns::Tools::to_value("2 % 3 : mondo", t);
     *
     * auto s = Patterns::Tools::to_string(t);
     * std::cout << s; // will print "2 % 3 : mondo""
     * @endcode
     *
     * See the documentation of the class Patterns::Tools::Convert, and of the
     * helper class Patterns::Tools::RankInfo for details on the separators you
     * should use in your string patterns when converting from a string to a
     * container type.
     *
     * Notice that the current content of variable @p t is ignored. Its type is
     * used to infer how to interpret the string. If the string is successfully
     * parsed, then @p t will be set to the parsed content of @p s.
     *
     * @author Luca Heltai, 2018
     */
    template <typename T>
    void
    to_value(const std::string &s, T &t);

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception.
     */
    DeclException2(ExcNoMatch,
                   std::string,
                   std::string,
                   << "The string " << arg1 << " does not match the pattern \""
                   << arg2 << "\"");
    //@}
  } // namespace Tools
} // namespace Patterns


// ---------------------- inline and template functions --------------------
namespace Patterns
{
  template <class... PatternTypes>
  Tuple::Tuple(const char *separator, const PatternTypes &... ps)
    : // forward to the version with std::string argument
    Tuple(std::string(separator), ps...)
  {}



  template <class... PatternTypes>
  Tuple::Tuple(const std::string &separator, const PatternTypes &... ps)
    : separator(separator)
  {
    static_assert(is_base_of_all<PatternBase, PatternTypes...>::value,
                  "Not all of the input arguments of this function "
                  "are derived from PatternBase");
    static_assert(sizeof...(ps) > 0,
                  "The number of PatternTypes must be greater than zero!");
    const auto pattern_pointers = {(static_cast<const PatternBase *>(&ps))...};
    for (const auto p : pattern_pointers)
      patterns.push_back(p->clone());
  }



  template <class... PatternTypes>
  Tuple::Tuple(const PatternTypes &... ps)
    : // forward to the version with the separator argument
    Tuple(std::string(":"), ps...)
  {}



  namespace Tools
  {
    namespace internal
    {
      /**
       * Store information about the rank types of the given class.
       *
       * A class has Rank equal to the number of different separators
       * that are required to uniquely identify its element(s) in a string.
       *
       * This class is used to detect whether the class T is compatible
       * with a Patterns::List pattern or with a Patterns::Map pattern.
       *
       * Objects like Point() or std::complex<double> are vector-likes, and
       * have vector_rank 1. Elementary types, like `int`, `unsigned int`,
       * `double`, etc. have vector_rank 0. `std::vector`, `std::list` and in
       * general containers have rank equal to 1 + vector_rank of the contained
       * type. Similarly for map types.
       *
       * A class with list_rank::value = 0 is either elementary or a
       * map. A class with map_rank::value = 0 is either a List compatible
       * class, or an elementary type.
       *
       * Elementary types are not compatible with Patterns::List, but non
       * elementary types, like Point(), or std::complex<double>, are compatible
       * with the List type. Adding more compatible types is a matter of adding
       * a specialization of this struct for the given type.
       *
       * @author Luca Heltai, 2017
       */
      template <class T, class Enable = void>
      struct RankInfo
      {
        static constexpr int list_rank = 0;
        static constexpr int map_rank  = 0;
      };
    } // namespace internal

    // Arithmetic types
    template <class T>
    struct Convert<T,
                   typename std::enable_if<std::is_arithmetic<T>::value>::type>
    {
      template <typename Dummy = T>
      static
        typename std::enable_if<std::is_same<Dummy, T>::value &&
                                  std::is_same<T, bool>::value,
                                std::unique_ptr<Patterns::PatternBase>>::type
        to_pattern()
      {
        return std_cxx14::make_unique<Patterns::Bool>();
      }

      template <typename Dummy = T>
      static
        typename std::enable_if<std::is_same<Dummy, T>::value &&
                                  !std::is_same<T, bool>::value &&
                                  std::is_integral<T>::value,
                                std::unique_ptr<Patterns::PatternBase>>::type
        to_pattern()
      {
        return std_cxx14::make_unique<Patterns::Integer>(
          std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max());
      }

      template <typename Dummy = T>
      static
        typename std::enable_if<std::is_same<Dummy, T>::value &&
                                  !std::is_same<T, bool>::value &&
                                  std::is_floating_point<T>::value,
                                std::unique_ptr<Patterns::PatternBase>>::type
        to_pattern()
      {
        return std_cxx14::make_unique<Patterns::Double>(
          std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max());
      }

      static std::string
      to_string(const T &                                     value,
                const std::unique_ptr<Patterns::PatternBase> &p =
                  Convert<T>::to_pattern())
      {
        std::stringstream str;
        if (std::is_same<T, unsigned char>::value ||
            std::is_same<T, signed char>::value || std::is_same<T, char>::value)
          str << static_cast<int>(value);
        else if (std::is_same<T, bool>::value)
          str << (static_cast<bool>(value) ? "true" : "false");
        else
          str << value;
        AssertThrow(p->match(str.str()),
                    ExcNoMatch(str.str(), p->description()));
        return str.str();
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &p =
                 Convert<T>::to_pattern())
      {
        AssertThrow(p->match(s), ExcNoMatch(s, p->description()));
        T value;
        if (std::is_same<T, bool>::value)
          value = (s == "true");
        else
          {
            std::istringstream is(s);
            if (std::is_same<T, unsigned char>::value ||
                std::is_same<T, signed char>::value ||
                std::is_same<T, char>::value)
              {
                int i;
                is >> i;
                value = i;
              }
            else
              is >> value;

            // If someone passes "123 abc" to the function, the method yields an
            // integer 123 alright, but the space terminates the read from the
            // string although there is more to come. This case, however, is
            // checked for in the call p->match(s) at the beginning of this
            // function, and would throw earlier. Here it is safe to assume that
            // if we didn't fail the conversion with the operator >>, then we
            // are good to go.
            AssertThrow(
              !is.fail(),
              ExcMessage("Failed to convert from \"" + s + "\" to the type \"" +
                         boost::core::demangle(typeid(T).name()) + "\""));
          }
        return value;
      }
    };

    namespace internal
    {
      constexpr std::array<const char *, 4> default_list_separator{
        {",", ";", "|", "%"}};
      constexpr std::array<const char *, 4> default_map_separator{
        {":", "=", "@", "#"}};

      // specialize a type for all of the STL containers and maps
      template <typename T>
      struct is_list_compatible : std::false_type
      {};
      template <typename T, std::size_t N>
      struct is_list_compatible<std::array<T, N>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::vector<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::deque<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::list<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::set<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::multiset<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::unordered_set<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::unordered_multiset<Args...>>
        : std::true_type
      {};

      template <typename T>
      struct is_map_compatible : std::false_type
      {};
      template <class Key, class T, class Compare, class Allocator>
      struct is_map_compatible<std::map<Key, T, Compare, Allocator>>
        : std::true_type
      {};
      template <class Key, class T, class Compare, class Allocator>
      struct is_map_compatible<std::multimap<Key, T, Compare, Allocator>>
        : std::true_type
      {};
      template <class Key, class T, class Hash, class KeyEqual, class Allocator>
      struct is_map_compatible<
        std::unordered_map<Key, T, Hash, KeyEqual, Allocator>> : std::true_type
      {};
      template <class Key, class T, class Hash, class KeyEqual, class Allocator>
      struct is_map_compatible<
        std::unordered_multimap<Key, T, Hash, KeyEqual, Allocator>>
        : std::true_type
      {};
    } // namespace internal

    // type trait to use the implementation type traits as well as decay the
    // type
    template <typename T>
    struct is_list_compatible
    {
      static constexpr bool const value =
        internal::is_list_compatible<typename std::decay<T>::type>::value;
    };

    template <typename T>
    struct is_map_compatible
    {
      static constexpr bool const value =
        internal::is_map_compatible<typename std::decay<T>::type>::value;
    };

    namespace internal
    {
      // Helper function for list_rank
      template <class T>
      constexpr int
      max_list_rank()
      {
        return RankInfo<T>::list_rank;
      }

      template <class T1, class T2, class... Types>
      constexpr int
      max_list_rank()
      {
        return std_cxx14::max(RankInfo<T1>::list_rank,
                              max_list_rank<T2, Types...>());
      }

      // Helper function for map_rank
      template <class T>
      constexpr int
      max_map_rank()
      {
        return RankInfo<T>::map_rank;
      }

      template <class T1, class T2, class... Types>
      constexpr int
      max_map_rank()
      {
        return std_cxx14::max(RankInfo<T1>::map_rank,
                              max_map_rank<T2, Types...>());
      }

      // Rank of vector types
      template <class T>
      struct RankInfo<
        T,
        typename std::enable_if<is_list_compatible<T>::value>::type>
      {
        static constexpr int list_rank =
          RankInfo<typename T::value_type>::list_rank + 1;
        static constexpr int map_rank =
          RankInfo<typename T::value_type>::map_rank;
      };

      // Rank of map types
      template <class T>
      struct RankInfo<
        T,
        typename std::enable_if<is_map_compatible<T>::value>::type>
      {
        static constexpr int list_rank =
          max_list_rank<typename T::key_type, typename T::mapped_type>() + 1;
        static constexpr int map_rank =
          max_map_rank<typename T::key_type, typename T::mapped_type>() + 1;
      };

      // Rank of Tensor types
      template <int rank, int dim, class Number>
      struct RankInfo<Tensor<rank, dim, Number>>
      {
        static constexpr int list_rank = rank + RankInfo<Number>::list_rank;
        static constexpr int map_rank  = RankInfo<Number>::map_rank;
      };

      template <int dim, class Number>
      struct RankInfo<Point<dim, Number>> : RankInfo<Tensor<1, dim, Number>>
      {};

      // Rank of complex types
      template <class Number>
      struct RankInfo<std::complex<Number>>
      {
        static constexpr int list_rank = RankInfo<Number>::list_rank + 1;
        static constexpr int map_rank  = RankInfo<Number>::map_rank;
      };

      // Rank of FunctionParser
      template <int dim>
      struct RankInfo<std::unique_ptr<FunctionParser<dim>>>
      {
        static constexpr int list_rank = 1;
        static constexpr int map_rank  = 0;
      };

      // Rank of ComponentMask
      template <>
      struct RankInfo<ComponentMask>
      {
        static constexpr int list_rank = 1;
        static constexpr int map_rank  = 0;
      };

      // Rank of std::pair
      template <class Key, class Value>
      struct RankInfo<std::pair<Key, Value>>
      {
        static constexpr int list_rank =
          std_cxx14::max(RankInfo<Key>::list_rank, RankInfo<Value>::list_rank);
        static constexpr int map_rank =
          std_cxx14::max(RankInfo<Key>::map_rank, RankInfo<Value>::map_rank) +
          1;
      };


      template <class... Types>
      struct RankInfo<std::tuple<Types...>>
      {
        static constexpr int list_rank = max_list_rank<Types...>();
        static constexpr int map_rank  = max_map_rank<Types...>() + 1;
      };
    } // namespace internal

    // stl containers
    template <class T>
    struct Convert<T,
                   typename std::enable_if<is_list_compatible<T>::value>::type>
    {
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");
        return std_cxx14::make_unique<Patterns::List>(
          *Convert<typename T::value_type>::to_pattern(),
          0,
          std::numeric_limits<unsigned int>::max(),
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1]);
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::List *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a "
                               "string to a List type."));
        auto                     base_p = p->get_base_pattern().clone();
        std::vector<std::string> vec(t.size());

        unsigned int i = 0;
        for (const auto &entry : t)
          vec[i++] = Convert<typename T::value_type>::to_string(entry, base_p);

        std::string s;
        if (vec.size() > 0)
          s = vec[0];
        for (unsigned int i = 1; i < vec.size(); ++i)
          s += p->get_separator() + " " + vec[i];

        AssertThrow(pattern->match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        AssertThrow(pattern->match(s), ExcNoMatch(s, pattern->description()));

        auto p = dynamic_cast<const Patterns::List *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List type."));

        auto base_p = p->get_base_pattern().clone();
        T    t;

        auto v = Utilities::split_string_list(s, p->get_separator());
        for (const auto &str : v)
          t.insert(t.end(),
                   Convert<typename T::value_type>::to_value(str, base_p));

        return t;
      }
    };

    // stl maps
    template <class T>
    struct Convert<T,
                   typename std::enable_if<is_map_compatible<T>::value>::type>
    {
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");
        static_assert(internal::RankInfo<T>::map_rank > 0,
                      "Cannot use this class for non Map-compatible types.");
        return std_cxx14::make_unique<Patterns::Map>(
          *Convert<typename T::key_type>::to_pattern(),
          *Convert<typename T::mapped_type>::to_pattern(),
          0,
          std::numeric_limits<unsigned int>::max(),
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1],
          internal::default_map_separator[internal::RankInfo<T>::map_rank - 1]);
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::Map *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a Map pattern to convert a string to "
                               "a Map compatible type."));
        auto                     key_p = p->get_key_pattern().clone();
        auto                     val_p = p->get_value_pattern().clone();
        std::vector<std::string> vec(t.size());

        unsigned int i = 0;
        for (const auto &ti : t)
          vec[i++] =
            Convert<typename T::key_type>::to_string(ti.first, key_p) +
            p->get_key_value_separator() +
            Convert<typename T::mapped_type>::to_string(ti.second, val_p);

        std::string s;
        if (vec.size() > 0)
          s = vec[0];
        for (unsigned int i = 1; i < vec.size(); ++i)
          s += p->get_separator() + " " + vec[i];

        AssertThrow(p->match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        AssertThrow(pattern->match(s), ExcNoMatch(s, pattern->description()));

        auto p = dynamic_cast<const Patterns::Map *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a Map pattern to convert a "
                               "string to a Map compatible type."));

        auto key_p = p->get_key_pattern().clone();
        auto val_p = p->get_value_pattern().clone();
        T    t;

        auto v = Utilities::split_string_list(s, p->get_separator());
        for (const auto &str : v)
          {
            auto key_val =
              Utilities::split_string_list(str, p->get_key_value_separator());
            AssertDimension(key_val.size(), 2);
            t.insert(std::make_pair(
              Convert<typename T::key_type>::to_value(key_val[0], key_p),
              Convert<typename T::mapped_type>::to_value(key_val[1])));
          }

        return t;
      }
    };

    // Tensors
    template <int rank, int dim, class Number>
    struct Convert<Tensor<rank, dim, Number>>
    {
      using T = Tensor<rank, dim, Number>;
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");
        return std_cxx14::make_unique<Patterns::List>(
          *Convert<typename T::value_type>::to_pattern(),
          dim,
          dim,
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1]);
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::List *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));
        auto                     base_p = p->get_base_pattern().clone();
        std::vector<std::string> vec(dim);

        for (unsigned int i = 0; i < dim; ++i)
          vec[i] = Convert<typename T::value_type>::to_string(t[i], base_p);

        std::string s;
        if (vec.size() > 0)
          s = vec[0];
        for (unsigned int i = 1; i < vec.size(); ++i)
          s += p->get_separator() + " " + vec[i];

        AssertThrow(p->match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        AssertThrow(pattern->match(s), ExcNoMatch(s, pattern->description()));

        auto p = dynamic_cast<const Patterns::List *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        auto base_p = p->get_base_pattern().clone();
        T    t;

        auto         v = Utilities::split_string_list(s, p->get_separator());
        unsigned int i = 0;
        for (const auto &str : v)
          t[i++] = Convert<typename T::value_type>::to_value(str, base_p);

        return t;
      }
    };

    // Points
    template <int dim, class Number>
    struct Convert<Point<dim, Number>>
    {
      using T = Point<dim, Number>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return Convert<Tensor<1, dim, Number>>::to_pattern();
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        return Convert<Tensor<1, dim, Number>>::to_string(
          Tensor<1, dim, Number>(t), pattern);
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        return T(Convert<Tensor<1, dim, Number>>::to_value(s, pattern));
      }
    };

    // Functions::FunctionParser
    template <int dim>
    struct Convert<std::unique_ptr<FunctionParser<dim>>>
    {
      using T = std::unique_ptr<FunctionParser<dim>>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");

        return std_cxx14::make_unique<Patterns::List>(
          Patterns::Anything(),
          1,
          Patterns::List::max_int_value,
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1]);
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::List *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        const auto &expressions = t->get_expressions();
        if (expressions.size() == 0)
          return std::string();

        std::string s = expressions[0];
        for (unsigned int i = 1; i < expressions.size(); ++i)
          s = s + p->get_separator() + expressions[i];

        AssertThrow(pattern->match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        AssertThrow(pattern->match(s), ExcNoMatch(s, pattern->description()));

        auto p = dynamic_cast<const Patterns::List *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        const auto expressions =
          Utilities::split_string_list(s, p->get_separator());

        T t = std_cxx14::make_unique<FunctionParser<dim>>(expressions.size());
        const std::string var =
          FunctionParser<dim>::default_variable_names() + ",t";
        const typename FunctionParser<dim>::ConstMap constants;
        t->initialize(var, expressions, constants, true);
        return t;
      }
    };

    // ComponentMask
    template <>
    struct Convert<ComponentMask>
    {
      using T = ComponentMask;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return Convert<std::vector<bool>>::to_pattern();
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        std::vector<bool> mask(t.size());
        for (unsigned int i = 0; i < t.size(); ++i)
          mask[i] = t[i];

        return Convert<std::vector<bool>>::to_string(mask, pattern);
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        const auto mask = Convert<std::vector<bool>>::to_value(s, pattern);
        return ComponentMask(mask);
      }
    };

    // Complex numbers
    template <class Number>
    struct Convert<std::complex<Number>>
    {
      using T = std::complex<Number>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");
        return std_cxx14::make_unique<Patterns::List>(
          *Convert<typename T::value_type>::to_pattern(),
          2,
          2,
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1]);
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::List *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        auto        base_p = p->get_base_pattern().clone();
        std::string s =
          Convert<typename T::value_type>::to_string(t.real(), base_p) +
          p->get_separator() + " " +
          Convert<typename T::value_type>::to_string(t.imag(), base_p);

        AssertThrow(pattern->match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      /**
       * Convert a string to a value, using the given pattern, or a default one.
       */
      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        AssertThrow(pattern->match(s), ExcNoMatch(s, pattern->description()));

        auto p = dynamic_cast<const Patterns::List *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        auto base_p = p->get_base_pattern().clone();

        auto v = Utilities::split_string_list(s, p->get_separator());
        AssertDimension(v.size(), 2);
        T t(Convert<typename T::value_type>::to_value(v[0], base_p),
            Convert<typename T::value_type>::to_value(v[1], base_p));
        return t;
      }
    };

    // Strings
    template <>
    struct Convert<std::string>
    {
      using T = std::string;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return std_cxx14::make_unique<Patterns::Anything>();
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        AssertThrow(pattern->match(t), ExcNoMatch(t, pattern->description()));
        return t;
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        AssertThrow(pattern->match(s), ExcNoMatch(s, pattern->description()));
        return s;
      }
    };

    // Pairs
    template <class Key, class Value>
    struct Convert<std::pair<Key, Value>>
    {
      using T = std::pair<Key, Value>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::map_rank > 0,
                      "Cannot use this class for non Map-compatible types.");
        return std_cxx14::make_unique<Patterns::Map>(
          *Convert<Key>::to_pattern(),
          *Convert<Value>::to_pattern(),
          1,
          1,
          // We keep the same list separator of the previous level, as this is
          // a map with only 1 possible entry
          internal::default_list_separator[internal::RankInfo<T>::list_rank],
          internal::default_map_separator[internal::RankInfo<T>::map_rank - 1]);
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        std::unordered_map<Key, Value> m;
        m.insert(t);
        std::string s = Convert<decltype(m)>::to_string(m, pattern);
        AssertThrow(pattern->match(s), ExcNoMatch(s, pattern->description()));
        return s;
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        std::unordered_map<Key, Value> m;
        m = Convert<decltype(m)>::to_value(s, pattern);
        return *m.begin();
      }
    };

    // Tuples
    template <class... Args>
    struct Convert<std::tuple<Args...>>
    {
      using T = std::tuple<Args...>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::map_rank > 0,
                      "Cannot use this class for non tuple-compatible types.");
        return std_cxx14::make_unique<Patterns::Tuple>(
          internal::default_map_separator[internal::RankInfo<T>::map_rank - 1],
          *Convert<Args>::to_pattern()...);
      }

      static std::string
      to_string(const T &                                     t,
                const std::unique_ptr<Patterns::PatternBase> &pattern =
                  Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::Tuple *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a Tuple pattern to convert a tuple "
                               "to a string."));

        const auto  string_array = Convert<T>::to_string_internal_2(t, *p);
        std::string str;
        for (unsigned int i = 0; i < string_array.size(); ++i)
          str += (i ? " " + p->get_separator() + " " : "") + string_array[i];
        AssertThrow(p->match(str), ExcNoMatch(str, p->description()));
        return str;
      }

      static T
      to_value(const std::string &                           s,
               const std::unique_ptr<Patterns::PatternBase> &pattern =
                 Convert<T>::to_pattern())
      {
        AssertThrow(pattern->match(s), ExcNoMatch(s, pattern->description()));

        auto p = dynamic_cast<const Patterns::Tuple *>(pattern.get());
        AssertThrow(p,
                    ExcMessage("I need a Tuple pattern to convert a string "
                               "to a tuple type."));

        auto v = Utilities::split_string_list(s, p->get_separator());

        return Convert<T>::to_value_internal_2(v, *p);
      }

    private:
      template <std::size_t... U>
      static std::array<std::string, std::tuple_size<T>::value>
      to_string_internal_1(const T &              t,
                           const Patterns::Tuple &pattern,
                           std_cxx14::index_sequence<U...>)
      {
        std::array<std::string, std::tuple_size<T>::value> a = {
          {Convert<typename std::tuple_element<U, T>::type>::to_string(
            std::get<U>(t), pattern.get_pattern(U).clone())...}};
        return a;
      }

      static std::array<std::string, std::tuple_size<T>::value>
      to_string_internal_2(const T &t, const Patterns::Tuple &pattern)
      {
        return Convert<T>::to_string_internal_1(
          t,
          pattern,
          std_cxx14::make_index_sequence<std::tuple_size<T>::value>{});
      }

      template <std::size_t... U>
      static T
      to_value_internal_1(const std::vector<std::string> &s,
                          const Patterns::Tuple &         pattern,
                          std_cxx14::index_sequence<U...>)
      {
        return std::make_tuple(
          Convert<typename std::tuple_element<U, T>::type>::to_value(
            s[U], pattern.get_pattern(U).clone())...);
      }

      static T
      to_value_internal_2(const std::vector<std::string> &s,
                          const Patterns::Tuple &         pattern)
      {
        return Convert<T>::to_value_internal_1(
          s,
          pattern,
          std_cxx14::make_index_sequence<std::tuple_size<T>::value>{});
      }
    };

    // Utility function with default Pattern
    template <typename T>
    std::string
    to_string(const T &t)
    {
      return Convert<T>::to_string(t);
    }

    // Utility function with default Pattern
    template <typename T>
    void
    to_value(const std::string &s, T &t)
    {
      t = Convert<T>::to_value(s);
    }
  } // namespace Tools
} // namespace Patterns


DEAL_II_NAMESPACE_CLOSE

#endif
