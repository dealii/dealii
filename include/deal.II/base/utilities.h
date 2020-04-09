// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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

#ifndef dealii_utilities_h
#define dealii_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <functional>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#ifdef DEAL_II_WITH_TRILINOS
#  include <Epetra_Comm.h>
#  include <Epetra_Map.h>
#  include <Teuchos_Comm.hpp>
#  include <Teuchos_RCP.hpp>
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif
#endif

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/core/demangle.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>

#ifdef DEAL_II_WITH_ZLIB
#  include <boost/iostreams/device/back_inserter.hpp>
#  include <boost/iostreams/filter/gzip.hpp>
#  include <boost/iostreams/filtering_stream.hpp>
#  include <boost/iostreams/stream.hpp>
#endif

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN

// forward declare Point
#ifndef DOXYGEN
template <int dim, typename Number>
class Point;
#endif

/**
 * A namespace for utility functions that are not particularly specific to
 * finite element computing or numerical programs, but nevertheless are needed
 * in various contexts when writing applications.
 *
 * @ingroup utilities
 * @author Wolfgang Bangerth, 2005
 */
namespace Utilities
{
  /**
   * Return a string of the form "deal.II version x.y.z" where "x.y.z"
   * identifies the version of deal.II you are using. This information
   * is also provided by the DEAL_II_PACKAGE_NAME and
   * DEAL_II_PACKAGE_VERSION preprocessor variables.
   */
  std::string
  dealii_version_string();

  /**
   * Assign to each point in @p points an index using the Hilbert space filling curve.
   * To that end, a bounding box for @p points will be determined, based on which their
   * integer coordinates are calculated.
   * The linear index is given as a dim-collection of bits, from high to low.
   * This is done in order to keep the maximum resolution in terms of bit depth
   * along each axis. Note that this dim-integer index can still be easily used
   * for sorting and ordering, for example using the lexicographic ordering of
   * tuples of integers.
   *
   * The depth of the Hilbert curve (i.e. the number of bits per dimension) by
   * default is equal to <code>64</code>.
   *
   * @note This function can also handle degenerate cases, e.g. when the bounding
   * box has zero size in one of the dimensions.
   */
  template <int dim, typename Number>
  std::vector<std::array<std::uint64_t, dim>>
  inverse_Hilbert_space_filling_curve(
    const std::vector<Point<dim, Number>> &points,
    const int                              bits_per_dim = 64);

  /**
   * Same as above, but for points in integer coordinates.
   */
  template <int dim>
  std::vector<std::array<std::uint64_t, dim>>
  inverse_Hilbert_space_filling_curve(
    const std::vector<std::array<std::uint64_t, dim>> &points,
    const int                                          bits_per_dim = 64);

  /**
   * Pack the least significant @p bits_per_dim bits from each element of @p index
   * (starting from last) into a single unsigned integer. The last element
   * of @p index will be used to set the first @p bits_per_dim bits in the
   * resulting integer, the second to last element is used to set the next @p bits_per_dim bits,
   * etc.. To fit all the data into the output, the following should hold
   * <code>bits\_per\_dim * dim <= 64</code>.
   *
   * The function is useful in debugging and visualization of indices returned
   * by inverse_Hilbert_space_filling_curve().
   *
   * @note There is no need to use this function in order to compare indices
   * returned by inverse_Hilbert_space_filling_curve(), as that can easily be
   * done via <code>std::lexicographical_compare()</code>.
   */
  template <int dim>
  std::uint64_t
  pack_integers(const std::array<std::uint64_t, dim> &index,
                const int                             bits_per_dim);

  /**
   * If the library is configured with ZLIB, then this function compresses the
   * input string and returns a non-zero terminated string containing the
   * compressed input.
   *
   * If the library was not configured with ZLIB enabled, the returned string
   * is identical to the input string.
   *
   * @param[in] input The string to compress
   *
   * @return A compressed version of the input string
   *
   * @authors Luca Heltai, Nicola Giuliani, 2020
   */
  std::string
  compress(const std::string &input);

  /**
   * If the library is configured with ZLIB, then this function assumes that the
   * input string has been compressed using the compress() function, and returns
   * the original decompressed string.
   *
   * If the library was not configured with ZLIB enabled, the returned string
   * is identical to the input string.
   *
   * @param[in] compressed_input A compressed string, as returned by the
   * function compress()
   *
   * @return The original uncompressed string.
   *
   * @authors Luca Heltai, Nicola Giuliani, 2020
   */
  std::string
  decompress(const std::string &compressed_input);

  /**
   * Encodes the binary input as a base64 string.
   *
   * Base64 is a group of binary-to-text encoding schemes that represent binary
   * data in an ASCII string format by translating it into a radix-64
   * representation. Base64 is designed to carry data stored in binary formats
   * across channels that only reliably support text content. It is used also
   * to store binary formats in a machine independent way.
   *
   * @param binary_input A vector of characters, representing your input as
   * binary data.
   * @return A string containing the binary input as a base64 string.
   *
   * @author Luca Heltai, 2020.
   */
  std::string
  encode_base64(const std::vector<unsigned char> &binary_input);

  /**
   * Decodes a base64 string into a binary output.
   *
   * This is the inverse of the encode_base64() function above.
   *
   * @param base64_input A string that contains the input in base64 format.
   * @return A vector of characters that represents your input as binary data.
   *
   * @author Luca Heltai, 2020.
   */
  std::vector<unsigned char>
  decode_base64(const std::string &base64_input);

  /**
   * Convert a number @p value to a string, with as many digits as given to
   * fill with leading zeros.
   *
   * If the second parameter is left at its default value, the number is not
   * padded with leading zeros. The result is then the same as if the standard
   * C++ `std::to_string` (or the older C function `itoa()`) had been called.
   *
   * This function takes an `unsigned int` as argument. As a consequence,
   * if you call it with a `signed int` (which is of course the same
   * type as `int`), the argument is implicitly converted to
   * unsigned integers and negative numbers may not be printed as you had
   * hoped. Similarly, if you call the function with a `long int`, the
   * printed result might show the effects of an overflow upon conversion
   * to `unsigned int`.
   *
   * @note The use of this function is discouraged and users should use
   * <code>Utilities::to_string()</code> instead. In its current
   * implementation the function simply calls <code>to_string@<unsigned
   * int@>()</code>.
   */
  std::string
  int_to_string(const unsigned int value,
                const unsigned int digits = numbers::invalid_unsigned_int);

  /**
   * Convert a number @p value to a string, with @p digits characters. The
   * string is padded with leading zeros, after a possible minus sign.
   * Therefore the total number of padding zeros is @p digits minus any signs,
   * decimal points and digits of @p value.
   *
   * If the second parameter is left at its default value, the number is not
   * padded with leading zeros. The result is then the same as if the C++
   * function `std::to_string()` had been called (for integral types),
   * or if `boost::lexical_cast()` had been called (for all other types).
   */
  template <typename number>
  std::string
  to_string(const number       value,
            const unsigned int digits = numbers::invalid_unsigned_int);

  /**
   * Determine how many digits are needed to represent numbers at most as
   * large as the given number.
   *
   * @author Niklas Fehn, 2019
   */
  unsigned int
  needed_digits(const unsigned int max_number);

  /**
   * This function allows to cut off a floating point number @p number
   * after @p n_digits of accuracy, i.e., after @p n_digits decimal places
   * in scientific floating point notation. When interpreted as rounding
   * operation, this function reduces the absolute value of a floating point
   * number and always rounds towards zero, since decimal places are simply
   * cut off.
   *
   * @author Niklas Fehn, 2019
   */
  template <typename Number>
  Number
  truncate_to_n_digits(const Number number, const unsigned int n_digits);

  /**
   * Given a string, convert it to an integer. Throw an assertion if that is
   * not possible.
   */
  int
  string_to_int(const std::string &s);

  /**
   * Return a string describing the dimensions of the object. Often, functions
   * in the deal.II library as well as in user codes need to define a string
   * containing the template dimensions of some objects defined using two
   * template parameters: dim (the topological dimension of the object) and
   * spacedim (the dimension of the embedding Euclidean space).  Since in all
   * deal.II classes, by default spacedim is equal to dimension, the above
   * string is usually contracted to "<dim>", instead of "<dim,spacedim>".
   * This function returns a string containing "dim" if dim is equal to
   * spacedim, otherwise it returns "dim,spacedim".
   */
  std::string
  dim_string(const int dim, const int spacedim);

  /**
   * Given a list of strings, convert it to a list of integers. Throw an
   * assertion if that is not possible.
   */
  std::vector<int>
  string_to_int(const std::vector<std::string> &s);

  /**
   * Given a string, convert it to an double. Throw an assertion if that is
   * not possible.
   */
  double
  string_to_double(const std::string &s);


  /**
   * Given a list of strings, convert it to a list of doubles. Throw an
   * assertion if that is not possible.
   */
  std::vector<double>
  string_to_double(const std::vector<std::string> &s);


  /**
   * Given a string that contains text separated by a @p delimiter, split it
   * into its components; for each component, remove leading and trailing
   * spaces. The default value of the delimiter is a comma, so that the
   * function splits comma separated lists of strings.
   *
   * To make data input from tables simpler, if the input string ends in a
   * delimiter (possibly followed by an arbitrary amount of whitespace), then
   * this last delimiter is ignored. For example,
   * @code
   *   Utilities::split_string_list("abc; def; ghi; ", ';');
   * @endcode
   * yields the same 3-element list of output <code>{"abc","def","ghi"}</code>
   * as you would get if the input had been
   * @code
   *   Utilities::split_string_list("abc; def; ghi", ';');
   * @endcode
   * or
   * @code
   *   Utilities::split_string_list("abc; def; ghi;", ';');
   * @endcode
   * As a consequence of this rule, a call like
   * @code
   *   Utilities::split_string_list(" ; ", ';');
   * @endcode
   * yields a one-element list. Because of the trimming of whitespace, the
   * single element is the empty string.
   *
   * This function can digest the case that the delimiter is a space. In this
   * case, it returns all words in the string. Combined with the rules above,
   * this implies that
   * @code
   *   Utilities::split_string_list("abc def ghi ", ' ');
   * @endcode
   * yields again the 3-element list of output
   * <code>{"abc","def","ghi"}</code> from above despite the presence of space
   * at the end of the string. Furthermore,
   * @code
   *   Utilities::split_string_list("      ", ' ');
   * @endcode
   * yields an empty list regardless of the number of spaces in the string.
   */
  std::vector<std::string>
  split_string_list(const std::string &s, const std::string &delimiter = ",");


  /**
   * Specialization of split_string_list() for the case where the delimiter
   * is a single char.
   */
  std::vector<std::string>
  split_string_list(const std::string &s, const char delimiter);


  /**
   * Take a text, usually a documentation or something, and try to break it
   * into individual lines of text at most @p width characters wide, by
   * breaking at positions marked by @p delimiter in the text. If this is not
   * possible, return the shortest lines that are longer than @p width.  The
   * default value of the delimiter is a space character. If original_text
   * contains newline characters (\n), the string is split at these locations,
   * too.
   */
  std::vector<std::string>
  break_text_into_lines(const std::string &original_text,
                        const unsigned int width,
                        const char         delimiter = ' ');

  /**
   * Return true if the given pattern string appears in the first position of
   * the string.
   */
  bool
  match_at_string_start(const std::string &name, const std::string &pattern);

  /**
   * Read a (signed) integer starting at the position in @p name indicated by
   * the second argument, and return this integer as a pair together with how
   * many characters it takes up in the string.
   *
   * If no integer can be read at the indicated position, return
   * (-1,numbers::invalid_unsigned_int)
   */
  std::pair<int, unsigned int>
  get_integer_at_position(const std::string &name, const unsigned int position);

  /**
   * Return a string with all occurrences of @p from in @p input replaced by
   * @p to.
   */
  std::string
  replace_in_string(const std::string &input,
                    const std::string &from,
                    const std::string &to);

  /**
   * Return a string with all standard whitespace characters (including
   * '<tt>\\t</tt>', '<tt>\\n</tt>', and '<tt>\\r</tt>') at the beginning and
   * end of @p input removed.
   */
  std::string
  trim(const std::string &input);

  /**
   * Generate a random number from a normalized Gaussian probability
   * distribution centered around @p a and with standard deviation @p sigma.
   * The returned number will be different every time the function is called.
   *
   * This function is reentrant, i.e., it can safely be called from multiple
   * threads at the same time. In addition, each thread will get the same
   * sequence of numbers every time. On the other hand, if you run
   * Threads::Task objects via the Threading Building Blocks, then tasks will
   * be assigned to mostly random threads, and may get a different sequence of
   * random numbers in different runs of the program, since a previous task
   * may already have consumed the first few random numbers generated for the
   * thread you're on. If this is a problem, you need to create your own
   * random number generator objects every time you want to start from a
   * defined point.
   *
   * @note Like the system function rand(), this function produces the same
   * sequence of random numbers every time a program is started. This is an
   * important property for debugging codes, but it makes it impossible to
   * really verify statistical properties of a code. For `rand()`, you can call
   * `srand()` to "seed" the random number generator to get different sequences
   * of random numbers every time a program is called. However, this function
   * does not allow seeding the random number generator. If you need this, as
   * above, use one of the C++ or BOOST facilities.
   */
  double
  generate_normal_random_number(const double a, const double sigma);

  /**
   * Return a string description of the type of the variable @p t.
   *
   * In general, C++ uses mangled names to identify types. This function
   * uses boost::core::demangle to return a human readable string describing
   * the type of the variable passed as argument.
   *
   * @author Luca Heltai, 2019.
   */
  template <class T>
  std::string
  type_to_string(const T &t);

  /**
   * Calculate a fixed power, provided as a template argument, of a number.
   *
   * This function provides an efficient way to calculate things like
   * <code>t^N</code> where <code>N</code> is a known number at compile time.
   *
   * Use this function as in <code>fixed_power@<dim@> (n)</code>.
   */
  template <int N, typename T>
  T
  fixed_power(const T t);

  /**
   * A replacement for <code>std::pow</code> that allows compile-time
   * calculations for constant expression arguments. The @p base must
   * be an integer type and the exponent @p iexp must not be negative.
   */
  template <typename T>
  constexpr T
  pow(const T base, const int iexp)
  {
#if defined(DEBUG) && defined(DEAL_II_HAVE_CXX14_CONSTEXPR)
    // Up to __builtin_expect this is the same code as in the 'Assert' macro.
    // The call to __builtin_expect turns out to be problematic.
    if (!(iexp >= 0))
      ::dealii::deal_II_exceptions::internals::issue_error_noreturn(
        ::dealii::deal_II_exceptions::internals::abort_or_throw_on_exception,
        __FILE__,
        __LINE__,
        __PRETTY_FUNCTION__,
        "iexp>=0",
        "ExcMessage(\"The exponent must not be negative!\")",
        ExcMessage("The exponent must not be negative!"));
#endif
    // The "exponentiation by squaring" algorithm used below has to be
    // compressed to one statement due to C++11's restrictions on constexpr
    // functions. A more descriptive version would be:
    //
    // <code>
    // if (iexp <= 0)
    //   return 1;
    //
    // // avoid overflow of one additional recursion with pow(base * base, 0)
    // if (iexp == 1)
    //   return base;
    //
    // // if the current exponent is not divisible by two,
    // // we need to account for that.
    // const unsigned int prefactor = (iexp % 2 == 1) ? base : 1;
    //
    // // a^b = (a*a)^(b/2)      for b even
    // // a^b = a*(a*a)^((b-1)/2 for b odd
    // return prefactor * dealii::Utilities::pow(base*base, iexp/2);
    // </code>

    static_assert(std::is_integral<T>::value, "Only integral types supported");

    return iexp <= 0 ?
             1 :
             (iexp == 1 ? base :
                          (((iexp % 2 == 1) ? base : 1) *
                           dealii::Utilities::pow(base * base, iexp / 2)));
  }

  /**
   * Optimized replacement for <tt>std::lower_bound</tt> for searching within
   * the range of column indices. Slashes execution time by approximately one
   * half for the present application, partly because the binary search is
   * replaced by a linear search for small loop lengths.
   *
   * Another reason for this function is rather obscure: when using the GCC
   * libstdc++ function std::lower_bound, complexity is O(log(N)) as required.
   * However, when using the debug version of the GCC libstdc++ as we do when
   * running the testsuite, then std::lower_bound tests whether the sequence
   * is in fact partitioned with respect to the pivot 'value' (i.e. in essence
   * that the sequence is sorted as required for binary search to work).
   * However, verifying this means that the complexity of std::lower_bound
   * jumps to O(N); we call this function O(N) times below, making the overall
   * complexity O(N**2). The consequence is that a few tests with big meshes
   * completely run off the wall time limit for tests and fail with the
   * libstdc++ debug mode
   *
   * This function simply makes the assumption that the sequence is sorted,
   * and we simply don't do the additional check.
   */
  template <typename Iterator, typename T>
  Iterator
  lower_bound(Iterator first, Iterator last, const T &val);


  /**
   * The same function as above, but taking an argument that is used to
   * compare individual elements of the sequence of objects pointed to by the
   * iterators.
   */
  template <typename Iterator, typename T, typename Comp>
  Iterator
  lower_bound(Iterator first, Iterator last, const T &val, const Comp comp);

  /**
   * Given a permutation vector (i.e. a vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for $i\neq j$), produce the reverse
   * permutation $q_i=N-1-p_i$.
   */
  template <typename Integer>
  std::vector<Integer>
  reverse_permutation(const std::vector<Integer> &permutation);

  /**
   * Given a permutation vector (i.e. a vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for $i\neq j$), produce the inverse
   * permutation $q_0\ldots q_{N-1}$ so that $q_{p_i}=p_{q_i}=i$.
   */
  template <typename Integer>
  std::vector<Integer>
  invert_permutation(const std::vector<Integer> &permutation);

  /**
   * Given an arbitrary object of type T, use boost::serialization utilities
   * to pack the object into a vector of characters and append it to the
   * given buffer. The number of elements that have been added to the buffer
   * will be returned. The object can be unpacked using the Utilities::unpack
   * function below.
   *
   * If the library has been compiled with ZLIB enabled, then the output buffer
   * can be compressed. This can be triggered with the parameter
   * @p allow_compression, and is only of effect if ZLIB is enabled.
   *
   * If many consecutive calls with the same buffer are considered, it is
   * recommended for reasons of performance to ensure that its capacity is
   * sufficient.
   *
   * @author Timo Heister, Wolfgang Bangerth, 2017.
   */
  template <typename T>
  size_t
  pack(const T &          object,
       std::vector<char> &dest_buffer,
       const bool         allow_compression = true);

  /**
   * Creates and returns a buffer solely for the given object, using the
   * above mentioned pack function.
   *
   * If the library has been compiled with ZLIB enabled, then the output buffer
   * can be compressed. This can be triggered with the parameter
   * @p allow_compression, and is only of effect if ZLIB is enabled.
   *
   * @author Timo Heister, Wolfgang Bangerth, 2017.
   */
  template <typename T>
  std::vector<char>
  pack(const T &object, const bool allow_compression = true);

  /**
   * Given a vector of characters, obtained through a call to the function
   * Utilities::pack, restore its content in an object of type T.
   *
   * This function uses boost::serialization utilities to unpack the object
   * from a vector of characters, and it is the inverse of the function
   * Utilities::pack().
   *
   * The @p allow_compression parameter denotes if the buffer to
   * read from could have been previously compressed with ZLIB, and
   * is only of effect if ZLIB is enabled.
   *
   * @note Since no arguments to this function depend on the template type
   *  @p T, you must manually specify the template argument when calling
   *  this function.
   *
   * @note If you want to pack() or unpack() arrays of objects, then the
   *  following works:
   *  @code
   *    double array[3] = {1,2,3};
   *    std::vector<char> buffer = Utilities::pack(array);
   *  @endcode
   *  However, the converse does not:
   *  @code
   *    array = Utilities::unpack<double[3]>(buffer);
   *  @endcode
   *  This is because C++ does not allow functions to return arrays.
   *  Consequently, there is a separate unpack() function for arrays, see
   *  below.
   *
   * @author Timo Heister, Wolfgang Bangerth, 2017.
   */
  template <typename T>
  T
  unpack(const std::vector<char> &buffer, const bool allow_compression = true);

  /**
   * Same unpack function as above, but takes constant iterators on
   * (a fraction of) a given packed buffer of type std::vector<char> instead.
   *
   * The @p allow_compression parameter denotes if the buffer to
   * read from could have been previously compressed with ZLIB, and
   * is only of effect if ZLIB is enabled.
   *
   * @author Timo Heister, Wolfgang Bangerth, 2017.
   */
  template <typename T>
  T
  unpack(const std::vector<char>::const_iterator &cbegin,
         const std::vector<char>::const_iterator &cend,
         const bool                               allow_compression = true);

  /**
   * Given a vector of characters, obtained through a call to the function
   * Utilities::pack, restore its content in an array of type T.
   *
   * This function uses boost::serialization utilities to unpack the object
   * from a vector of characters, and it is the inverse of the function
   * Utilities::pack().
   *
   * The @p allow_compression parameter denotes if the buffer to
   * read from could have been previously compressed with ZLIB, and
   * is only of effect if ZLIB is enabled.
   *
   * @note This function exists due to a quirk of C++. Specifically,
   *  if you want to pack() or unpack() arrays of objects, then the
   *  following works:
   *  @code
   *    double array[3] = {1,2,3};
   *    std::vector<char> buffer = Utilities::pack(array);
   *  @endcode
   *  However, the converse does not:
   *  @code
   *    array = Utilities::unpack<double[3]>(buffer);
   *  @endcode
   *  This is because C++ does not allow functions to return arrays.
   *  The current function therefore allows to write
   *  @code
   *    Utilities::unpack(buffer, array);
   *  @endcode
   *  Note that unlike the other unpack() function, it is not necessary
   *  to explicitly specify the template arguments since they can be
   *  deduced from the second argument.
   *
   * @author Timo Heister, Wolfgang Bangerth, 2017.
   */
  template <typename T, int N>
  void
  unpack(const std::vector<char> &buffer,
         T (&unpacked_object)[N],
         const bool allow_compression = true);

  /**
   * Same unpack function as above, but takes constant iterators on
   * (a fraction of) a given packed buffer of type std::vector<char> instead.
   *
   * The @p allow_compression parameter denotes if the buffer to
   * read from could have been previously compressed with ZLIB, and
   * is only of effect if ZLIB is enabled.
   *
   * @author Timo Heister, Wolfgang Bangerth, 2017.
   */
  template <typename T, int N>
  void
  unpack(const std::vector<char>::const_iterator &cbegin,
         const std::vector<char>::const_iterator &cend,
         T (&unpacked_object)[N],
         const bool allow_compression = true);

  /**
   * Convert an object of type `std::unique_ptr<From>` to an object of
   * type `std::unique_ptr<To>`, where it is assumed that we can cast
   * the pointer to `From` to a pointer to `To` using a `dynamic_cast`
   * -- in other words, we assume that `From` and `To` are connected
   * through a class hierarchy, and that the object pointed to is in
   * fact of a type that contains both a `From` and a `To`. An example
   * is if either `To` is derived from `From` or the other way around.
   *
   * The function throws an exception of type `std::bad_cast` if the
   * `dynamic_cast` does not succeed. This is the same exception you
   * would get if a regular `dynamic_cast` between object types (but not
   * pointer types) does not succeed.
   *
   * An example of how this function works is as follows:
   * @code
   *   // A base class. Assume that it has virtual
   *   // functions so that dynamic_cast can work.
   *   class B
   *   {
   *     ...
   *   };
   *
   *   // A derived class
   *   class D : public B
   *   {
   *     ...
   *   };
   *
   *   // A factory function
   *   std::unique_ptr<B> create_object (...)
   *   {
   *     ...
   *   }
   *
   *   void foo (...)
   *   {
   *     std::unique_ptr<B> b = create_object (...);
   *
   *     // Assume that we know for some reason that the object above must
   *     // have created a D object but returned it as a std::unique_ptr<B>.
   *     // In order to access the D functionality, we need to cast the
   *     // pointer. Use the equivalent to dynamic_cast:
   *     std::unique_ptr<D> d = dynamic_unique_cast<D>(std::move(b));
   *
   *     // If the object really was a D, then 'd' now points to it. Note
   *     // also that in accordance with the semantics of std::unique_ptr,
   *     // it was necessary to std::move the 'b' object, and indeed 'b'
   *     // now no longer points to anything -- ownership has been
   *     // transferred to 'd'!
   * @endcode
   *
   * @note This function does not try to convert the `Deleter` objects stored
   *   by `std::unique_ptr` objects. The function therefore only works if the
   *   deleter objects are at their defaults, i.e., if they are of type
   *   `std::default_delete<To>` and `std::default_delete<From>`.
   */
  template <typename To, typename From>
  std::unique_ptr<To>
  dynamic_unique_cast(std::unique_ptr<From> &&p)
  {
    // Let's see if we can cast from 'From' to 'To'. If so, do the cast,
    // and then release the pointer from the old
    // owner
    if (To *cast = dynamic_cast<To *>(p.get()))
      {
        std::unique_ptr<To> result(cast);
        p.release();
        return result;
      }
    else
      throw std::bad_cast();
  }

  /**
   * A namespace for utility functions that probe system properties.
   *
   * @ingroup utilities
   */
  namespace System
  {
    /**
     * Return the CPU load as returned by "uptime". Note that the
     * interpretation of this number depends on the actual number of
     * processors in the machine. This is presently only implemented on Linux,
     * using the /proc/loadavg pseudo-file, on other systems we simply return
     * zero.
     */
    double
    get_cpu_load();

    /**
     * Return the instruction set extension for vectorization as described by
     * DEAL_II_VECTORIZATION_WIDTH_IN_BITS in vectorization.h as a string. The
     * list of possible return values is:
     *
     * <table>
     * <tr>
     *   <td><tt>VECTORIZATION_LEVEL</tt></td>
     *   <td>Return Value</td>
     *   <td>Width in bits</td>
     * </tr>
     * <tr>
     *   <td>0</td>
     *   <td>disabled</td>
     *   <td>64</td>
     * </tr>
     * <tr>
     *   <td>1</td>
     *   <td>SSE2/AltiVec</td>
     *   <td>128</td>
     * </tr>
     * <tr>
     *   <td>2</td>
     *   <td>AVX</td>
     *   <td>256</td>
     * </tr>
     * <tr>
     *   <td>3</td>
     *   <td>AVX512</td>
     *   <td>512</td>
     * </tr>
     * </table>
     */
    const std::string
    get_current_vectorization_level();

    /**
     * Structure that hold information about memory usage in kB. Used by
     * get_memory_stats(). See man 5 proc entry /status for details.
     */
    struct MemoryStats
    {
      /**
       * Peak virtual memory size in kB.
       */
      unsigned long int VmPeak;

      /**
       * Current virtual memory size in kB.
       */
      unsigned long int VmSize;

      /**
       * Peak resident memory size in kB. Also known as "high water mark" (HWM).
       */
      unsigned long int VmHWM;

      /**
       * Current resident memory size in kB. Also known as "resident set size"
       * (RSS).
       */
      unsigned long int VmRSS;
    };


    /**
     * Fill the @p stats structure with information about the memory
     * consumption of this process. This is only implemented on Linux.
     */
    void
    get_memory_stats(MemoryStats &stats);


    /**
     * Return the name of the host this process runs on.
     */
    std::string
    get_hostname();


    /**
     * Return the present time as HH:MM:SS.
     */
    std::string
    get_time();

    /**
     * Return the present date as YYYY/MM/DD. MM and DD may be either one or
     * two digits.
     */
    std::string
    get_date();

    /**
     * Call the system function posix_memalign, or a replacement function if
     * not available, to allocate memory with a certain minimal alignment. The
     * first argument will then return a pointer to this memory block that can
     * be released later on through a standard <code>free</code> call.
     *
     * @param memptr The address of a pointer variable that will after this
     * call point to the allocated memory.
     * @param alignment The minimal alignment of the memory block, in bytes.
     * @param size The size of the memory block to be allocated, in bytes.
     *
     * @note This function checks internally for error codes, rather than
     * leaving this task to the calling site.
     */
    void
    posix_memalign(void **memptr, std::size_t alignment, std::size_t size);
  } // namespace System


#ifdef DEAL_II_WITH_TRILINOS
  /**
   * This namespace provides some of the basic structures used in the
   * initialization of the Trilinos objects (e.g., matrices, vectors, and
   * preconditioners).
   */
  namespace Trilinos
  {
    /**
     * Return a Trilinos Epetra_Comm object needed for creation of
     * Epetra_Maps.
     *
     * If deal.II has been configured to use a compiler that does not support
     * MPI then the resulting communicator will be a serial one. Otherwise,
     * the communicator will correspond to MPI_COMM_WORLD, i.e. a communicator
     * that encompasses all processes within this MPI universe.
     */
    const Epetra_Comm &
    comm_world();

    /**
     * Return a Trilinos Epetra_Comm object needed for creation of
     * Epetra_Maps.
     *
     * If deal.II has been configured to use a compiler that does not support
     * MPI then the resulting communicator will be a serial one. Otherwise,
     * the communicator will correspond to MPI_COMM_SELF, i.e. a communicator
     * that comprises only this one processor.
     */
    const Epetra_Comm &
    comm_self();

    /**
     * Return a Teuchos::Comm object needed for creation of Tpetra::Maps.
     *
     * If deal.II has been configured to use a compiler that does not support
     * MPI then the resulting communicator will be a serial one. Otherwise,
     * the communicator will correspond to MPI_COMM_SELF, i.e. a communicator
     * that comprises only this one processor.
     */
    const Teuchos::RCP<const Teuchos::Comm<int>> &
    tpetra_comm_self();

    /**
     * Given a communicator, duplicate it. If the given communicator is
     * serial, that means to just return a copy of itself. On the other hand,
     * if it is %parallel, we duplicate the underlying MPI_Comm object: we
     * create a separate MPI communicator that contains the same processors
     * and in the same order but has a separate identifier distinct from the
     * given communicator. The function returns a pointer to a new object of a
     * class derived from Epetra_Comm. The caller of this function needs to
     * assume ownership of this function. The returned object should be
     * destroyed using the destroy_communicator() function.
     *
     * This facility is used to separate streams of communication. For
     * example, a program could simply use MPI_Comm_World for everything. But
     * it is easy to come up with scenarios where sometimes not all processors
     * participate in a communication that is intended to be global -- for
     * example if we assemble a matrix on a coarse mesh with fewer cells than
     * there are processors, some processors may not sync their matrices with
     * the rest because they haven't written into it because they own no
     * cells. That's clearly a bug. However, if these processors just continue
     * their work, and the next %parallel operation happens to be a sync on a
     * different matrix, then the sync could succeed -- by accident, since
     * different processors are talking about different matrices.
     *
     * This kind of situation can be avoided if we use different communicators
     * for different matrices which reduces the likelihood that communications
     * meant to be separate aren't recognized as such just because they happen
     * on the same communicator. In addition, it is conceivable that some MPI
     * operations can be parallelized using multiple threads because their
     * communicators identifies the communication in question, not their
     * relative timing as is the case in a sequential program that just uses a
     * single communicator.
     */
    Epetra_Comm *
    duplicate_communicator(const Epetra_Comm &communicator);

    /**
     * Given an Epetra communicator that was created by the
     * duplicate_communicator() function, destroy the underlying MPI
     * communicator object and reset the Epetra_Comm object to a the result of
     * comm_self().
     *
     * It is necessary to call this function at the time when the result of
     * duplicate_communicator() is no longer needed. The reason is that in
     * that function, we first create a new MPI_Comm object and then create an
     * Epetra_Comm around it. While we can take care of destroying the latter,
     * it doesn't destroy the communicator since it can only assume that it
     * may also be still used by other objects in the program. Consequently,
     * we have to take care of destroying it ourselves, explicitly.
     *
     * This function does exactly that. Because this has to happen while the
     * Epetra_Comm object is still around, it first resets the latter and then
     * destroys the communicator object.
     *
     * @note If you call this function on an Epetra_Comm object that is not
     * created by duplicate_communicator(), you are likely doing something
     * quite wrong. Don't do this.
     */
    void
    destroy_communicator(Epetra_Comm &communicator);

    /**
     * Return the number of MPI processes there exist in the given
     * @ref GlossMPICommunicator "communicator"
     * object. If this is a sequential job (i.e., the program
     * is not using MPI at all, or is using MPI but has been started with
     * only one MPI process), then the communicator necessarily involves
     * only one process and the function returns 1.
     */
    unsigned int
    get_n_mpi_processes(const Epetra_Comm &mpi_communicator);

    /**
     * Return the number of the present MPI process in the space of processes
     * described by the given communicator. This will be a unique value for
     * each process between zero and (less than) the number of all processes
     * (given by get_n_mpi_processes()).
     */
    unsigned int
    get_this_mpi_process(const Epetra_Comm &mpi_communicator);

    /**
     * Given a Trilinos Epetra map, create a new map that has the same
     * subdivision of elements to processors but uses the given communicator
     * object instead of the one stored in the first argument. In essence,
     * this means that we create a map that communicates among the same
     * processors in the same way, but using a separate channel.
     *
     * This function is typically used with a communicator that has been
     * obtained by the duplicate_communicator() function.
     */
    Epetra_Map
    duplicate_map(const Epetra_BlockMap &map, const Epetra_Comm &comm);
  } // namespace Trilinos

#endif


} // namespace Utilities


// --------------------- inline functions

namespace Utilities
{
  template <int N, typename T>
  inline T
  fixed_power(const T x)
  {
    Assert(
      !std::is_integral<T>::value || (N >= 0),
      ExcMessage(
        "The non-type template parameter N must be a non-negative integer for integral type T"));

    if (N == 0)
      return T(1.);
    else if (N < 0)
      return T(1.) / fixed_power<-N>(x);
    else
      // Use exponentiation by squaring:
      return ((N % 2 == 1) ? x * fixed_power<N / 2>(x * x) :
                             fixed_power<N / 2>(x * x));
  }



  template <class T>
  inline std::string
  type_to_string(const T &t)
  {
    return boost::core::demangle(typeid(t).name());
  }



  template <typename Iterator, typename T>
  inline Iterator
  lower_bound(Iterator first, Iterator last, const T &val)
  {
    return Utilities::lower_bound(first, last, val, std::less<T>());
  }



  template <typename Iterator, typename T, typename Comp>
  inline Iterator
  lower_bound(Iterator first, Iterator last, const T &val, const Comp comp)
  {
    // verify that the two iterators are properly ordered. since
    // we need operator- for the iterator type anyway, do the
    // test as follows, rather than via 'last >= first'
    Assert(last - first >= 0,
           ExcMessage(
             "The given iterators do not satisfy the proper ordering."));

    unsigned int len = static_cast<unsigned int>(last - first);

    if (len == 0)
      return first;

    while (true)
      {
        // if length equals 8 or less,
        // then do a rolled out
        // search. use a switch without
        // breaks for that and roll-out
        // the loop somehow
        if (len < 8)
          {
            switch (len)
              {
                case 7:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 6:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 5:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 4:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 3:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 2:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 1:
                  if (!comp(*first, val))
                    return first;
                  return first + 1;
                default:
                  // indices seem
                  // to not be
                  // sorted
                  // correctly!? or
                  // did len
                  // become==0
                  // somehow? that
                  // shouldn't have
                  // happened
                  Assert(false, ExcInternalError());
              }
          }



        const unsigned int half   = len >> 1;
        const Iterator     middle = first + half;

        // if the value is larger than
        // that pointed to by the
        // middle pointer, then the
        // insertion point must be
        // right of it
        if (comp(*middle, val))
          {
            first = middle + 1;
            len -= half + 1;
          }
        else
          len = half;
      }
  }


  // --------------------- non-inline functions

  template <typename T>
  size_t
  pack(const T &          object,
       std::vector<char> &dest_buffer,
       const bool         allow_compression)
  {
    // the data is never compressed when we can't use zlib.
    (void)allow_compression;

    std::size_t size = 0;

    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
    //
    // we have to work around the fact that GCC 4.8.x claims to be C++
    // conforming, but is not actually as it does not implement
    // std::is_trivially_copyable.
#if __GNUG__ && __GNUC__ < 5
    if (__has_trivial_copy(T) && sizeof(T) < 256)
#else
#  ifdef DEAL_II_WITH_CXX17
    if constexpr (std::is_trivially_copyable<T>() && sizeof(T) < 256)
#  else
    if (std::is_trivially_copyable<T>() && sizeof(T) < 256)
#  endif
#endif
      {
        const std::size_t previous_size = dest_buffer.size();
        dest_buffer.resize(previous_size + sizeof(T));

        std::memcpy(dest_buffer.data() + previous_size, &object, sizeof(T));

        size = sizeof(T);
      }
    else
      {
        // use buffer as the target of a compressing
        // stream into which we serialize the current object
        const std::size_t previous_size = dest_buffer.size();
#ifdef DEAL_II_WITH_ZLIB
        if (allow_compression)
          {
            boost::iostreams::filtering_ostream out;
            out.push(
              boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(
                boost::iostreams::gzip::default_compression)));
            out.push(boost::iostreams::back_inserter(dest_buffer));

            boost::archive::binary_oarchive archive(out);
            archive << object;
            out.flush();
          }
        else
#endif
          {
            std::ostringstream              out;
            boost::archive::binary_oarchive archive(out);
            archive << object;

            const std::string &s = out.str();
            dest_buffer.reserve(dest_buffer.size() + s.size());
            std::move(s.begin(), s.end(), std::back_inserter(dest_buffer));
          }

        size = dest_buffer.size() - previous_size;
      }

    return size;
  }


  template <typename T>
  std::vector<char>
  pack(const T &object, const bool allow_compression)
  {
    std::vector<char> buffer;
    pack<T>(object, buffer, allow_compression);
    return buffer;
  }


  template <typename T>
  T
  unpack(const std::vector<char>::const_iterator &cbegin,
         const std::vector<char>::const_iterator &cend,
         const bool                               allow_compression)
  {
    T object;

    // the data is never compressed when we can't use zlib.
    (void)allow_compression;

    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
    //
    // we have to work around the fact that GCC 4.8.x claims to be C++
    // conforming, but is not actually as it does not implement
    // std::is_trivially_copyable.
#if __GNUG__ && __GNUC__ < 5
    if (__has_trivial_copy(T) && sizeof(T) < 256)
#else
#  ifdef DEAL_II_WITH_CXX17
    if constexpr (std::is_trivially_copyable<T>() && sizeof(T) < 256)
#  else
    if (std::is_trivially_copyable<T>() && sizeof(T) < 256)
#  endif
#endif
      {
        Assert(std::distance(cbegin, cend) == sizeof(T), ExcInternalError());
        std::memcpy(&object, &*cbegin, sizeof(T));
      }
    else
      {
        std::string decompressed_buffer;

        // first decompress the buffer
#ifdef DEAL_II_WITH_ZLIB
        if (allow_compression)
          {
            boost::iostreams::filtering_ostream decompressing_stream;
            decompressing_stream.push(boost::iostreams::gzip_decompressor());
            decompressing_stream.push(
              boost::iostreams::back_inserter(decompressed_buffer));
            decompressing_stream.write(&*cbegin, std::distance(cbegin, cend));
          }
        else
#endif
          {
            decompressed_buffer.assign(cbegin, cend);
          }

        // then restore the object from the buffer
        std::istringstream              in(decompressed_buffer);
        boost::archive::binary_iarchive archive(in);

        archive >> object;
      }

    return object;
  }


  template <typename T>
  T
  unpack(const std::vector<char> &buffer, const bool allow_compression)
  {
    return unpack<T>(buffer.cbegin(), buffer.cend(), allow_compression);
  }


  template <typename T, int N>
  void
  unpack(const std::vector<char>::const_iterator &cbegin,
         const std::vector<char>::const_iterator &cend,
         T (&unpacked_object)[N],
         const bool allow_compression)
  {
    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
    //
    // we have to work around the fact that GCC 4.8.x claims to be C++
    // conforming, but is not actually as it does not implement
    // std::is_trivially_copyable.
    if (
#if __GNUG__ && __GNUC__ < 5
      __has_trivial_copy(T)
#else
      std::is_trivially_copyable<T>()
#endif
      && sizeof(T) * N < 256)
      {
        Assert(std::distance(cbegin, cend) == sizeof(T) * N,
               ExcInternalError());
        std::memcpy(unpacked_object, &*cbegin, sizeof(T) * N);
      }
    else
      {
        std::string decompressed_buffer;

        // first decompress the buffer
        (void)allow_compression;
#ifdef DEAL_II_WITH_ZLIB
        if (allow_compression)
          {
            boost::iostreams::filtering_ostream decompressing_stream;
            decompressing_stream.push(boost::iostreams::gzip_decompressor());
            decompressing_stream.push(
              boost::iostreams::back_inserter(decompressed_buffer));
            decompressing_stream.write(&*cbegin, std::distance(cbegin, cend));
          }
        else
#endif
          {
            decompressed_buffer.assign(cbegin, cend);
          }

        // then restore the object from the buffer
        std::istringstream              in(decompressed_buffer);
        boost::archive::binary_iarchive archive(in);

        archive >> unpacked_object;
      }
  }


  template <typename T, int N>
  void
  unpack(const std::vector<char> &buffer,
         T (&unpacked_object)[N],
         const bool allow_compression)
  {
    unpack<T, N>(buffer.cbegin(),
                 buffer.cend(),
                 unpacked_object,
                 allow_compression);
  }



  template <typename Integer>
  std::vector<Integer>
  reverse_permutation(const std::vector<Integer> &permutation)
  {
    const std::size_t n = permutation.size();

    std::vector<Integer> out(n);
    for (std::size_t i = 0; i < n; ++i)
      out[i] = n - 1 - permutation[i];

    return out;
  }



  template <typename Integer>
  std::vector<Integer>
  invert_permutation(const std::vector<Integer> &permutation)
  {
    const std::size_t n = permutation.size();

    std::vector<Integer> out(n, numbers::invalid_unsigned_int);

    for (std::size_t i = 0; i < n; ++i)
      {
        AssertIndexRange(permutation[i], n);
        out[permutation[i]] = i;
      }

    // check that we have actually reached
    // all indices
    for (std::size_t i = 0; i < n; ++i)
      Assert(out[i] != numbers::invalid_unsigned_int,
             ExcMessage("The given input permutation had duplicate entries!"));

    return out;
  }
} // namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#ifndef DOXYGEN
namespace boost
{
  namespace serialization
  {
    // Provides boost and c++11 with a way to serialize tuples and pairs
    // automatically
    template <int N>
    struct Serialize
    {
      template <class Archive, typename... Args>
      static void
      serialize(Archive &ar, std::tuple<Args...> &t, const unsigned int version)
      {
        ar &std::get<N - 1>(t);
        Serialize<N - 1>::serialize(ar, t, version);
      }
    };

    template <>
    struct Serialize<0>
    {
      template <class Archive, typename... Args>
      static void
      serialize(Archive &ar, std::tuple<Args...> &t, const unsigned int version)
      {
        (void)ar;
        (void)t;
        (void)version;
      }
    };

    template <class Archive, typename... Args>
    void
    serialize(Archive &ar, std::tuple<Args...> &t, const unsigned int version)
    {
      Serialize<sizeof...(Args)>::serialize(ar, t, version);
    }
  } // namespace serialization
} // namespace boost
#endif

#endif
