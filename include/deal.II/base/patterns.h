// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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

#ifndef dealii__patterns_h
#define dealii__patterns_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/archive/basic_archive.hpp>
#include <boost/property_tree/ptree_serialization.hpp>

#include <map>
#include <vector>
#include <string>
#include <memory>

DEAL_II_NAMESPACE_OPEN

// forward declarations for interfaces and friendship
class LogStream;
class MultipleParameterLoop;

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
    virtual ~PatternBase ();

    /**
     * Return <tt>true</tt> if the given string matches the pattern.
     */
    virtual bool match (const std::string &test_string) const = 0;

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
    virtual std::string description (const OutputStyle style=Machine) const = 0;

    /**
     * Return a pointer to an exact copy of the object. This is necessary
     * since we want to store objects of this type in containers, were we need
     * to copy objects without knowledge of their actual data type (we only
     * have pointers to the base class).
     *
     * Ownership of the objects returned by this function is passed to the
     * caller of this function.
     */
    virtual std::unique_ptr<PatternBase> clone () const = 0;

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
    virtual std::size_t memory_consumption () const;
  };

  /**
   * Return pointer to the correct derived class based on description.
   */
  std::unique_ptr<PatternBase> pattern_factory (const std::string &description);

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
     */
    Integer (const int lower_bound = min_int_value,
             const int upper_bound = max_int_value);

    /**
     * Return <tt>true</tt> if the string is an integer and its value is
     * within the specified range.
     */
    virtual bool match (const std::string &test_string) const;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. If bounds were specified to the constructor, then include them
     * into this description.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * Creates new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Integer> create (const std::string &description);

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
    Double (const double lower_bound = min_double_value,
            const double upper_bound = max_double_value);

    /**
     * Return <tt>true</tt> if the string is a number and its value is within
     * the specified range.
     */
    virtual bool match (const std::string &test_string) const;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. If bounds were specified to the constructor, then include them
     * into this description.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * Creates a new object on the heap using @p new if the given
     * @p description is a valid format (for example created by calling
     * description() on an existing object), or @p nullptr otherwise. Ownership
     * of the returned object is transferred to the caller of this function,
     * which should be freed using @p delete.
     */
    static std::unique_ptr<Double> create(const std::string &description);

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
    Selection (const std::string &seq);

    /**
     * Return <tt>true</tt> if the string is an element of the description
     * list passed to the constructor.
     */
    virtual bool match (const std::string &test_string) const;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the list of valid strings passed to the
     * constructor.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t memory_consumption () const;

    /**
     * Creates new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Selection> create (const std::string &description);

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
    List (const PatternBase  &base_pattern,
          const unsigned int  min_elements = 0,
          const unsigned int  max_elements = max_int_value,
          const std::string  &separator = ",");


    /**
     * Return the internally stored separator.
     */
    const std::string &get_separator() const;

    /**
     * Return the internally stored base pattern.
     */
    const PatternBase &get_base_pattern() const;

    /**
     * Copy constructor.
     */
    List (const List &other);

    /**
     * Return <tt>true</tt> if the string is a comma-separated list of strings
     * each of which match the pattern given to the constructor.
     */
    virtual bool match (const std::string &test_string) const;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * Creates new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<List> create (const std::string &description);

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t memory_consumption () const;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception.
     */
    DeclException2 (ExcInvalidRange,
                    int, int,
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
   * between pairs other than the comma, and a delimeter between key and value
   * other than column.
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
    Map (const PatternBase  &key_pattern,
         const PatternBase  &value_pattern,
         const unsigned int  min_elements = 0,
         const unsigned int  max_elements = max_int_value,
         const std::string  &separator = ",",
         const std::string  &key_value_separator = ":");

    /**
     * Copy constructor.
     */
    Map (const Map &other);

    /**
     * Return <tt>true</tt> if the string is a comma-separated list of strings
     * each of which match the pattern given to the constructor.
     */
    virtual bool match (const std::string &test_string) const;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * Creates new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Map> create (const std::string &description);

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t memory_consumption () const;

    /**
     * Return a reference to the key pattern.
     */
    const PatternBase &get_key_pattern() const;

    /**
     * Return a reference to the value pattern.
     */
    const PatternBase &get_value_pattern() const;

    /**
     * Return the separator of the map entries.
     */
    const std::string &get_separator() const;

    /**
     * Return the key-value separator.
     */
    const std::string &get_key_value_separator() const;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception.
     */
    DeclException2 (ExcInvalidRange,
                    int, int,
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
    MultipleSelection (const std::string &seq);

    /**
     * Return <tt>true</tt> if the string is an element of the description
     * list passed to the constructor.
     */
    virtual bool match (const std::string &test_string) const;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the list of valid strings passed to the
     * constructor.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * Creates new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<MultipleSelection> create (const std::string &description);

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t memory_consumption () const;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception.
     */
    DeclException1 (ExcCommasNotAllowed,
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
    Bool ();

    /**
     * Return a description of the pattern that valid strings are expected to
     * match.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * Creates new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Bool> create(const std::string &description);

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
    Anything ();

    /**
     * Return <tt>true</tt> if the string matches its constraints, i.e.
     * always.
     */
    virtual bool match (const std::string &test_string) const;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the string <tt>"[Anything]"</tt>.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * Creates new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<Anything> create(const std::string &description);

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
    FileName (const FileType type = input);

    /**
     * Return <tt>true</tt> if the string matches its constraints, i.e.
     * always.
     */
    virtual bool match (const std::string &test_string) const;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the string <tt>"[Filename]"</tt>.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * file type flag
     */
    FileType  file_type;

    /**
     * Creates new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<FileName> create (const std::string &description);

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
    DirectoryName ();

    /**
     * Return <tt>true</tt> if the string matches its constraints, i.e.
     * always.
     */
    virtual bool match (const std::string &test_string) const;

    /**
     * Return a description of the pattern that valid strings are expected to
     * match. Here, this is the string <tt>"[Filename]"</tt>.
     */
    virtual std::string description (const OutputStyle style=Machine) const;

    /**
     * Return a copy of the present object, which is newly allocated on the
     * heap. Ownership of that object is transferred to the caller of this
     * function.
     */
    virtual std::unique_ptr<PatternBase> clone () const;

    /**
     * Creates new object if the start of description matches
     * description_init.  Ownership of that object is transferred to the
     * caller of this function.
     */
    static std::unique_ptr<DirectoryName> create(const std::string &description);

  private:
    /**
     * Initial part of description
     */
    static const char *description_init;
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif
