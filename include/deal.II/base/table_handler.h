// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_table_handler_h
#define dealii_table_handler_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <boost/serialization/map.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <cstdint>
#include <map>
#include <ostream>
#include <string>
#include <variant>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class TableHandler;
#endif

namespace internal
{
  /**
   * A <tt>TableEntry</tt> stores the value of a table entry. It can either be
   * of type int, unsigned int, std::uint64_t, double or std::string. In
   * essence, this structure is the same as `std::variant<int,unsigned
   * int,std::uint64_t,double,std::string>` but we wrap this object in a
   * structure for which we can write a function that can serialize it. This is
   * also why the function is not in fact of type std::any.
   */
  struct TableEntry
  {
  public:
    /**
     * Default constructor.
     */
    TableEntry() = default;

    /**
     * Constructor. Initialize this table element with the value
     * <code>t</code>.
     */
    template <typename T>
    explicit TableEntry(const T &t);

    /**
     * Return the value stored by this object. The template type T must be one
     * of <code>int,unsigned int,std::uint64_t,double,std::string</code> and it
     * must match the data type of the object originally stored in this
     * TableEntry object.
     */
    template <typename T>
    T
    get() const;

    /**
     * Return the numeric value of this object if data has been stored in it
     * either as an integer, an unsigned integer,std::uint64_t, or a double.
     *
     * @return double
     */
    double
    get_numeric_value() const;

    /**
     * Cache the contained value with the given formatting and return it. The
     * given parameters from the column definition are used for the
     * formatting. The value is cached as a string internally in cached_value.
     * The cache needs to be invalidated with this routine if the formatting
     * of the column changes.
     */
    void
    cache_string(bool scientific, unsigned int precision) const;

    /**
     * Return the value cached using cache_string(). This is just a wrapper
     * around cached_value.
     */
    const std::string &
    get_cached_string() const;


    /**
     * Return a TableEntry object that has the same data type of the stored
     * value but with a value that is default constructed for this data type.
     * This is used to pad columns below previously set ones.
     */
    TableEntry
    get_default_constructed_copy() const;

    /**
     * Write the data of this object to a stream for the purpose of
     * serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    save(Archive &ar, const unsigned int version) const;

    /**
     * Read the data of this object from a stream for the purpose of
     * serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
    /**
     * Write and read the data of this object from a stream for the purpose
     * of serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    serialize(Archive &archive, const unsigned int version);
#else
    // This macro defines the serialize() method that is compatible with
    // the templated save() and load() method that have been implemented.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  private:
    /**
     * Abbreviation for the data type stored by this object.
     */
    using value_type =
      std::variant<int, unsigned int, std::uint64_t, double, std::string>;

    /**
     * Stored value.
     */
    value_type value;

    /**
     * Cache the current value as a string.
     */
    mutable std::string cached_value;

    friend class dealii::TableHandler;
  };
} // namespace internal


/**
 * The TableHandler stores TableEntries of arbitrary value type and writes the
 * table as text or in tex format to an output stream. The value type actually
 * may vary from column to column and from row to row.
 *
 * <h3>Usage</h3>
 *
 * The most important function is the templatized function
 * <code>add_value(const std::string &key, const T value)</code> that adds a
 * column with the name <tt>key</tt> to the table if this column does not yet
 * exist and adds the given value of type <tt>T</tt> (which must be one of
 * <tt>int</tt>, <tt>unsigned int</tt>, <tt>double</tt>, <tt>std::string</tt>)
 * to this column.  After the table is complete there are different
 * possibilities of output, e.g., into a latex file with write_tex() or as
 * text with write_text().
 *
 * Two (or more) columns may be merged into a "supercolumn" by twice (or
 * multiple) calling add_column_to_supercolumn(), see there. Additionally
 * there is a function to set for each column the precision of the output of
 * numbers, and there are several functions to prescribe the format and the
 * captions the columns are written with in tex mode.
 *
 * A detailed explanation of this class is also given in the step-13 tutorial
 * program.
 *
 *
 * <h3>Example</h3>
 *
 * This is a simple example demonstrating the usage of this class. The first
 * column includes the numbers $i=1 \dots n$, the second $1^2 \dots n^2$, the
 * third $\sqrt{1}\dots\sqrt{n}$, where the second and third columns are merged
 * into one supercolumn with the superkey <tt>squares and roots</tt>.
 * Additionally the first column is aligned to the right (the default was
 * <tt>centered</tt>) and the precision of the square roots are set to be 6
 * (instead of 4 as default).
 *
 * @code
 * TableHandler table;
 * for (unsigned int i = 1; i <= n; ++i)
 *   {
 *     table.add_value("numbers", i);
 *     table.add_value("squares", i * i);
 *     table.add_value("square roots", std::sqrt(i));
 *   }
 * // merge the second and third column
 * table.add_column_to_supercolumn("squares", "squares and roots");
 * table.add_column_to_supercolumn("square roots", "squares and roots");
 *
 * // additional settings
 * table.set_tex_format("numbers", "r");
 * table.set_precision("square roots", 6);
 *
 * // output
 * std::ofstream out_file("number_table.tex");
 * table.write_tex(out_file);
 * out_file.close();
 * @endcode
 *
 *
 * <h3>Dealing with sparse data: auto-fill mode</h3>
 *
 * When generating output, TableHandler expects that all columns have the
 * exact same number of elements in it so that the result is in fact a table.
 * This assumes that in each of the iterations (time steps, nonlinear
 * iterations, etc) you fill every single column. On the other hand, this may
 * not always be what you want to do. For example, it could be that the
 * function that computes the nonlinear residual is only called every few time
 * steps; or, a function computing statistics of the mesh is only called
 * whenever the mesh is in fact refined. In these cases, the add_value()
 * function will be called less often for some columns and the column would
 * therefore have fewer elements; furthermore, these elements would not be
 * aligned with the rows that contain the other data elements that were
 * produced during this iteration. An entirely different scenario is that the
 * table is filled and at a later time we use the data in there to compute the
 * elements of other rows; the ConvergenceTable class does something like
 * this.
 *
 * To support both scenarios, the TableHandler class has a property called <i
 * >auto-fill mode</i>. By default, auto-fill mode is off, but it can be
 * enabled by calling set_auto_fill_mode(). If auto-fill mode is enabled we
 * use the following algorithm:
 *
 * - When calling <code>add_value(key, value)</code>, we count the number of
 * elements in the column corresponding to <code>key</code>. Let's call this
 * number $m$.
 *
 * - We also determine the maximal number of elements in the other columns;
 * call it $n$.
 *
 * - If $m < n-1$ then we add $n-m-1$ copies of the object <code>T()</code> to
 * this column. Here, <code>T</code> is the data type of the given
 * <code>value</code>. For example, if <code>T</code> is a numeric type, then
 * <code>T()</code> is the number zero; if <code>T</code> is
 * <code>std::string</code>, then <code>T()</code> is the empty string
 * <code>""</code>.
 *
 * - Add the given value to this column.
 *
 * Padding the column with default elements makes sure that after the addition
 * the column has as many entries as the longest other column. In other words,
 * if we have skipped previous invocations of add_value() for a given key,
 * then the padding will enter default values into this column.
 *
 * The algorithm as described will fail if you try to skip adding values for a
 * key if adding an element for this key is the first thing you want to do for
 * a given iteration or time step, since we would then pad to the length of
 * the longest column of the <i>previous</i> iteration or time step. You may
 * have to re-order adding to this column to a different spot in your program,
 * after adding to a column that will always be added to; or, you may want to
 * start every iteration by adding the number of the iteration to the table,
 * for example in column 1.
 *
 * In the case above, we have always padded columns <b>above</b> the element
 * that is being added to a column. However, there is also a case where we
 * have to pad <b>below</b>. Namely, if a previous row has been completely
 * filled using TableHandler::add_value(), subsequent rows have been filled
 * partially, and we then ask for output via write_text() or write_tex(). In
 * that case, the last few rows that have been filled only partially need to
 * be padded below the last element that has been added to them. As before, we
 * do that by using default constructed objects of the same type as the last
 * element of that column.
 *
 * @ingroup textoutput
 */
class TableHandler
{
public:
  /**
   * Set of options how a table should be formatted when output with the
   * write_text() function. The following possibilities exist:
   *
   * - <code>table_with_headers</code>: The table is formatted in such a way
   * that the contents are aligned under the key of each column, i.e. the key
   * sits atop each column. This is suitable for tables with few columns where
   * the entire table can be displayed on the screen. Output looks like this:
   *   @code
   *     key1 key2 key3
   *     0    0    ""
   *     1    0    ""
   *     2    13   a
   *     1    0    ""
   *   @endcode
   * - <code>table_with_separate_column_description</code>: This is a better
   * format when there are many columns and the table as a whole can not be
   * displayed on the screen. Here, the column keys are first listed
   * one-by-one on lines of their own, and are numbered for better readability.
   * In addition, each of these description lines are prefixed by '#' to mark
   * these lines as comments for programs that want to read the following
   * table as data and should ignore these descriptive lines. GNUPLOT is one
   * such program that will automatically ignore lines so prefixed. Output
   * with this option looks like this:
   *   @code
   *     # 1: key1
   *     # 2: key2
   *     # 3: key3
   *     0 0  ""
   *     1 0  ""
   *     2 13 a
   *     1 0  ""
   *   @endcode
   * - <code>simple_table_with_separate_column_description</code>: This format
   * is very similar to <code>table_with_separate_column_description</code>,
   * but it skips aligning the columns with additional white space. This
   * increases the performance of write_text() for large tables. Example
   * output:
   *   @code
   *     # 1: key1
   *     # 2: key2
   *     # 3: key3
   *     0 0 ""
   *     1 0 ""
   *     2 13 a
   *     1 0 ""
   *   @endcode
   * - <code>org_mode_table</code>: Outputs to org-mode (http://orgmode.org/)
   * table format. It is easy to convert org-mode tables to HTML/LaTeX/csv.
   * Example output:
   *   @code
   *   | key1 | key2 | key3 |
   *   | 0    | 0    | ""   |
   *   | 1    | 0    | ""   |
   *   | 2    | 13   | a    |
   *   | 1    | 0    | ""   |
   *   @endcode
   */
  enum TextOutputFormat
  {
    /**
     * Print the table with headers.
     */
    table_with_headers,
    /**
     * Print the table with separate lines for each column label.
     */
    table_with_separate_column_description,
    /**
     * Like table_with_separate_column_description, but without aligning the
     * column containing the column labels.
     */
    simple_table_with_separate_column_description,
    /**
     * Print the table in org mode format.
     */
    org_mode_table
  };

  /**
   * Constructor.
   */
  TableHandler();


  /**
   * Declare the existence of a column in the table by giving it a name.
   * As discussed in the documentation of the class, this is not usually
   * necessary -- just adding a value for a given column key via the
   * add_value() function also declares the column. This function is
   * therefore only necessary in cases where you want a column to
   * also show up even if you never add an entry to any row in this column;
   * or, more likely, if you want to prescribe the order in which columns
   * are later printed by declaring columns in a particular order before
   * entries are ever put into them.
   *
   * (The latter objective can also be achieved by adding entries to
   * the table in whatever order they are produced by a program,
   * and later calling set_column_order(). However, this approach
   * requires knowing -- in one central place of your software --
   * all of the columns keys that other parts of the software have
   * written into, and how they should be sorted. This is easily
   * possible for small programs, but may not be feasible for
   * large code bases in which parts of the code base are only
   * executed based on run-time parameters.)
   */
  void
  declare_column(const std::string &key);

  /**
   * Adds a column (if not yet existent) with the key <tt>key</tt> and adds
   * the value of type <tt>T</tt> to the column. Values of type <tt>T</tt>
   * must be convertible to one of <code>int, unsigned int, double,
   * std::uint64_t, std::string</code> or a compiler error will result.
   */
  template <typename T>
  void
  add_value(const std::string &key, const T value);

  /**
   * If a row is only partially filled, then set all elements of that
   * row for which no elements exist in a particular column to the
   * empty string. This is akin to the 'auto_fill_mode' described in
   * the introduction, but more general because it allows you to start
   * writing into a column for a new row without having to know that
   * that column had been written to in the previous row.
   *
   * If all columns have been written into in the current row, then
   * this function doesn't do anything at all. In other words,
   * conceptually the function "completes" the current row, though its
   * use case is to "start" a new row.
   */
  void
  start_new_row();

  /**
   * Switch auto-fill mode on or off. See the general documentation of this
   * class for a description of what auto-fill mode does.
   */
  void
  set_auto_fill_mode(const bool state);

  /**
   * Creates a supercolumn (if not yet existent) and includes column to it.
   * The keys of the column and the supercolumn are <tt>key</tt> and
   * <tt>superkey</tt>, respectively.  To merge two columns <tt>c1</tt> and
   * <tt>c2</tt> to a supercolumn <tt>sc</tt> hence call
   * <tt>add_column_to_supercolumn(c1,sc)</tt> and
   * <tt>add_column_to_supercolumn(c2,sc)</tt>.
   *
   * Concerning the order of the columns, the supercolumn replaces the first
   * column that is added to the supercolumn. Within the supercolumn the order
   * of output follows the order the columns are added to the supercolumn.
   */
  void
  add_column_to_supercolumn(const std::string &key,
                            const std::string &superkey);

  /**
   * Change the order of columns and supercolumns in the table.
   *
   * <tt>new_order</tt> includes the keys and superkeys of the columns and
   * supercolumns in the order the user would like them to be output. If a
   * superkey is included the keys of the subcolumns need not be explicitly
   * mentioned in this vector.  The order of subcolumns within a supercolumn
   * is not changeable and remains in the order in which the columns are added
   * to the supercolumn.
   *
   * This function may also be used to break big tables with too many columns
   * into smaller ones. For example, you can call this function with the first
   * five columns and then call one of the <tt>write_*</tt> functions, then
   * call this function with the next five columns and again <tt>write_*</tt>,
   * and so on.
   */
  void
  set_column_order(const std::vector<std::string> &new_order);

  /**
   * Set the <tt>precision</tt> e.g. double or float variables are written
   * with. <tt>precision</tt> is the same as in calling
   * <tt>out<<setprecision(precision)</tt>.
   */
  void
  set_precision(const std::string &key, const unsigned int precision);

  /**
   * Set the <tt>scientific_flag</tt>. True means scientific, false means
   * fixed point notation.
   */
  void
  set_scientific(const std::string &key, const bool scientific);

  /**
   * Set the caption of the column <tt>key</tt> for tex output. You may want
   * to chose this different from <tt>key</tt>, if it contains formulas or
   * similar constructs.
   */
  void
  set_tex_caption(const std::string &key, const std::string &tex_caption);

  /**
   * Set the tex caption of the entire <tt>table</tt> for tex output.
   */
  void
  set_tex_table_caption(const std::string &table_caption);

  /**
   * Set the label of this <tt>table</tt> for tex output.
   */
  void
  set_tex_table_label(const std::string &table_label);

  /**
   * Set the caption the supercolumn <tt>superkey</tt> for tex output.
   * You may want to chose this different from <tt>superkey</tt>, if it
   * contains formulas or similar constructs.
   */
  void
  set_tex_supercaption(const std::string &superkey,
                       const std::string &tex_supercaption);

  /**
   * Set the tex output format of a column, e.g. <tt>c</tt>, <tt>r</tt>,
   * <tt>l</tt>, or <tt>p{3cm}</tt>. The default is <tt>c</tt>. Also if this
   * function is not called for a column, the default is preset to be
   * <tt>c</tt>.
   */
  void
  set_tex_format(const std::string &key, const std::string &format = "c");

  /**
   * Write table as formatted text to the given stream. The text is formatted
   * in such as way that it represents data as formatted columns of text. To
   * avoid problems when reading these tables automatically, for example for
   * postprocessing, if an entry in a cell of this table is empty (i.e. it has
   * been created by calling the add_value() function with an empty string),
   * then the entry of the table is printed as <code>""</code>.
   *
   * The second argument indicates how column keys are to be displayed. See
   * the description of TextOutputFormat for more information.
   */
  void
  write_text(std::ostream          &out,
             const TextOutputFormat format = table_with_headers) const;

  /**
   * Write table as a tex file. If @p with_header is set to false, then no
   * <code>\\documentclass{...}</code>, <code>\\begin{document}</code> and
   * <code>\\end{document}</code> are used. In this way the file can be
   * included into an existing tex file using a command like
   * <code>\\input{table_file}</code>.
   */
  void
  write_tex(std::ostream &file, const bool with_header = true) const;

  /**
   * Clear the rows of the table, i.e. calls clear() on all the underlying
   * storage data structures.
   */
  void
  clear();

  /**
   * Remove all values added at the current row. This is useful when, for
   * example, a time-step is rejected and all data recorded about it needs to
   * be discarded.
   */
  void
  clear_current_row();

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception
   */
  DeclException1(ExcColumnNotExistent,
                 std::string,
                 << "Column <" << arg1 << "> does not exist.");

  /**
   * Exception
   */
  DeclException1(ExcSuperColumnNotExistent,
                 std::string,
                 << "Supercolumn <" << arg1 << "> does not exist.");

  /**
   * Exception
   */
  DeclException1(ExcColumnOrSuperColumnNotExistent,
                 std::string,
                 << "Column or supercolumn <" << arg1 << "> does not exist.");

  /**
   * Exception
   */
  DeclException4(ExcWrongNumberOfDataEntries,
                 std::string,
                 int,
                 std::string,
                 int,
                 << "Column <" << arg1 << "> has " << arg2
                 << " rows, but Column <" << arg3 << "> has " << arg4
                 << " rows.");

  /**
   * Exception
   */
  DeclException1(ExcUndefinedTexFormat,
                 std::string,
                 << '<' << arg1 << "> is not a tex column format. Use "
                 << "'l', 'c', or 'r' to indicate left, centered, or "
                 << "right aligned text.");
  /** @} */
protected:
  /**
   * Structure encapsulating all the data that is needed to describe one
   * column of a table.
   */
  struct Column
  {
    /**
     * Constructor needed by <tt>std::map</tt>.
     */
    Column();

    /**
     * Constructor.
     */
    explicit Column(const std::string &tex_caption);

    /**
     * Pad this column with default constructed elements to the number of rows
     * given by the argument.
     */
    void
    pad_column_below(const unsigned int length);

    /**
     * Write the data of this object to a stream for the purpose of
     * serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    save(Archive &ar, const unsigned int version) const;

    /**
     * Read the data of this object from a stream for the purpose of
     * serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
    /**
     * Write and read the data of this object from a stream for the purpose
     * of serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
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
     * Invalidates the string cache of all the entries and recomputes the
     * maximum length max_length.
     */
    void
    invalidate_cache();

    /**
     * List of entries within this column. Values are always immediately
     * converted to strings to provide a uniform method of lookup.
     */
    std::vector<internal::TableEntry> entries;

    /**
     * The caption of the column in tex output.  By default, this is the key
     * string that is given to the <tt>TableHandler</tt> by
     * <tt>TableHandler::add_value(...)</tt>. This may be changed by calling
     * <tt>TableHandler::set_tex_caption(...)</tt>.
     */
    std::string tex_caption;

    /**
     * The column format in tex output.  By default, this is <tt>"c"</tt>,
     * meaning `centered'. This may be changed by calling
     * <tt>TableHandler::set_tex_format(...)</tt> with <tt>"c", "r", "l"</tt>
     * for centered, right or left.
     */

    std::string tex_format;

    /**
     * Double or float entries are written with this precision (set by the
     * user).  The default is 4.
     */
    unsigned int precision;

    /**
     * <tt>scientific</tt>=false means fixed point notation.
     */
    bool scientific;

    /**
     * Flag that may be used by derived classes for arbitrary purposes.
     *
     * In particular, the ConvergenceTable class uses the flag to denote
     * columns for which convergence information has already been computed, or
     * should not be computed at all.
     */
    unsigned int flag;

    /**
     * This entry caches the maximum length in characters for all entries in
     * this table.
     */
    unsigned int max_length;
  };

  /**
   * Help function that gives a vector of the keys of all columns that are
   * mentioned in <tt>column_order</tt>, where each supercolumn key is
   * replaced by its subcolumn keys.
   *
   * This function implicitly checks the consistency of the data. The result
   * is returned in <tt>sel_columns</tt>.
   */
  void
  get_selected_columns(std::vector<std::string> &sel_columns) const;

  /**
   * Builtin function, that gives the number of rows in the table and that
   * checks if the number of rows is equal in every column. This function is
   * e.g. called before writing output.
   */
  unsigned int
  n_rows() const;

  /**
   * A variable storing the column and supercolumn keys in the order desired by
   * the user. By default this is the order of adding the columns. This order
   * may be changed by set_column_order().
   */
  std::vector<std::string> column_order;

  /**
   * A map from the column keys to the columns (not supercolumns).
   *
   * The field is declared mutable so that the write_text() and write_tex()
   * functions can be const, even though they may pad columns below if
   * 'auto_fill_mode' is on.
   */
  mutable std::map<std::string, Column> columns;

  /**
   * A map from each supercolumn key to the keys of its subcolumns in the right
   * order.  It is allowed that a supercolumn has got the same key as a
   * column.
   *
   * Note that we do not use a <tt>multimap</tt> here since the order of
   * column keys for each supercolumn key is relevant.
   */
  std::map<std::string, std::vector<std::string>> supercolumns;

  /**
   * A map from the supercolumn keys to the captions of the supercolumns that
   * are used in tex output.
   *
   * By default these are just the supercolumn keys but they may be changed by
   * <tt>set_tex_supercaptions(...)</tt>.
   */
  std::map<std::string, std::string> tex_supercaptions;

  /**
   * The caption of the table itself.
   */
  std::string tex_table_caption;
  /**
   * The label of the table.
   */
  std::string tex_table_label;

  /**
   * Flag indicating whether auto-fill mode should be used.
   */
  bool auto_fill_mode;
};


namespace internal
{
  template <typename T>
  TableEntry::TableEntry(const T &t)
    : value(t)
  {}


  template <typename T>
  T
  TableEntry::get() const
  {
    // we don't quite know the data type in 'value', but
    // it must be one of the ones in the type list of the
    // std::variant. so if T is not in the list, or if
    // the data stored in the TableEntry is not of type
    // T, then we will get an exception that we can
    // catch and produce an error message
    try
      {
        return std::get<T>(value);
      }
    catch (...)
      {
        Assert(false,
               ExcMessage(
                 "This TableEntry object does not store a datum of type T"));
        throw;
      }
  }



  template <class Archive>
  void
  TableEntry::save(Archive &ar, const unsigned int) const
  {
    // write first an identifier for the kind
    // of data stored and then the actual
    // data, in its correct data type
    if (std::holds_alternative<int>(value))
      {
        const int p = std::get<int>(value);
        char      c = 'i';
        ar &c    &p;
      }
    else if (std::holds_alternative<unsigned int>(value))
      {
        const unsigned int p = std::get<unsigned int>(value);
        char               c = 'u';
        ar &c             &p;
      }
    else if (std::holds_alternative<double>(value))
      {
        const double p = std::get<double>(value);
        char         c = 'd';
        ar &c       &p;
      }
    else if (std::holds_alternative<std::string>(value))
      {
        const std::string p = std::get<std::string>(value);
        char              c = 's';
        ar &c            &p;
      }
    else if (std::holds_alternative<std::uint64_t>(value))
      {
        const std::uint64_t p = std::get<std::uint64_t>(value);
        char                c = 'l';
        ar &c              &p;
      }
    else
      DEAL_II_ASSERT_UNREACHABLE();
  }



  template <class Archive>
  void
  TableEntry::load(Archive &ar, const unsigned int)
  {
    // following what we do in the save()
    // function, first read in the data type
    // as a one-character id, and then read
    // the data
    char c;
    ar  &c;

    switch (c)
      {
        case 'i':
          {
            int val;
            ar &val;
            value = val;
            break;
          }

        case 'u':
          {
            unsigned int val;
            ar          &val;
            value = val;
            break;
          }

        case 'd':
          {
            double val;
            ar    &val;
            value = val;
            break;
          }

        case 's':
          {
            std::string val;
            ar         &val;
            value = val;
            break;
          }

        case 'l':
          {
            std::uint64_t val;
            ar           &val;
            value = val;
            break;
          }

        default:
          DEAL_II_ASSERT_UNREACHABLE();
      }
  }
} // namespace internal



template <typename T>
void
TableHandler::add_value(const std::string &key, const T value)
{
  // see if the column already exists
  if (columns.find(key) == columns.end())
    declare_column(key);

  if (auto_fill_mode == true)
    {
      // follow the algorithm given in the introduction to this class
      // of padding columns as necessary
      unsigned int max_col_length = 0;
      for (const auto &column : columns)
        max_col_length =
          std::max(max_col_length,
                   static_cast<unsigned int>(column.second.entries.size()));

      while (columns[key].entries.size() + 1 < max_col_length)
        {
          columns[key].entries.push_back(internal::TableEntry(T()));
          const internal::TableEntry &entry = columns[key].entries.back();
          entry.cache_string(columns[key].scientific, columns[key].precision);
          columns[key].max_length =
            std::max(columns[key].max_length,
                     static_cast<unsigned int>(
                       entry.get_cached_string().size()));
        }
    }

  // now push the value given to this function
  columns[key].entries.push_back(internal::TableEntry(value));
  const internal::TableEntry &entry = columns[key].entries.back();
  entry.cache_string(columns[key].scientific, columns[key].precision);
  columns[key].max_length =
    std::max(columns[key].max_length,
             static_cast<unsigned int>(entry.get_cached_string().size()));
}



template <class Archive>
void
TableHandler::Column::save(Archive &ar, const unsigned int /*version*/) const
{
  ar &entries &tex_caption &tex_format &precision &scientific &flag &max_length;
}



template <class Archive>
void
TableHandler::Column::load(Archive &ar, const unsigned int /*version*/)
{
  ar &entries &tex_caption &tex_format &precision &scientific &flag &max_length;
  invalidate_cache();
}


template <class Archive>
void
TableHandler::serialize(Archive &ar, const unsigned int)
{
  ar &column_order &columns &supercolumns &tex_supercaptions &tex_table_caption
    &tex_table_label &auto_fill_mode;
}


DEAL_II_NAMESPACE_CLOSE

#endif
