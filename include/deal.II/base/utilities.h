// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_utilities_h
#define dealii_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/types.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/core/demangle.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#ifdef DEAL_II_WITH_ZLIB
#  include <boost/iostreams/filter/gzip.hpp>
#endif

DEAL_II_NAMESPACE_OPEN

// forward declare Point
#ifndef DOXYGEN
template <int dim, typename Number>
DEAL_II_CXX20_REQUIRES(dim >= 0)
class Point;
#endif

/**
 * A namespace for utility functions that are not particularly specific to
 * finite element computing or numerical programs, but nevertheless are needed
 * in various contexts when writing applications.
 *
 * @ingroup utilities
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
   * While the function takes the argument `t`, it does not actually use
   * its value but only the type of `t` for its output.
   */
  template <class T>
  std::string
  type_to_string(const T &t);

  /**
   * Calculate a fixed power, provided as a template argument, of a number.
   *
   * This function provides an efficient way to calculate things like
   * $t^N$ where <code>N</code> is a known number at compile time.
   * The function computes the power of $t$ via the "recursive doubling"
   * approach in which, for example, $t^7$ is computed as
   * @f{align*}{
   *   t^7 = (tttt)(tt)(t)
   * @f}
   * where computing $tt$ requires one product, computing $tttt$ is
   * achieved by multiplying the previously computed $tt$ by itself
   * (requiring another multiplication), and then the product is computed
   * via two more multiplications for a total of 4 multiplications
   * instead of the naively necessary 6.
   *
   * The major savings this function generates result, however, from the
   * fact that it exploits that we have an integer power of the argument $t$.
   * The alternative to computing such powers, `std::pow(t,7)` uses the
   * `std::pow` function that takes the exponent as a floating point number
   * and, because it has to cater to the complexities of the general situation,
   * is vastly slower.
   *
   * Use this function as in `fixed_power<dim> (t)` or `fixed_power<7> (t)`.
   */
  template <int N, typename T>
  constexpr T
  fixed_power(const T t);

  /**
   * A replacement for <code>std::pow</code> that allows compile-time
   * calculations for constant expression arguments and if the exponent
   * is a positive integer. The @p base must be an arithmetic type
   * (i.e., an integer or floating point type),
   * and the exponent @p iexp must not be negative.
   */
  template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
  constexpr DEAL_II_HOST_DEVICE T
  pow(const T base, const int iexp);

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
   * Given an arbitrary object of type `T`, use boost::serialization utilities
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
   * This function considers a number of special cases for which packing (and
   * unpacking) can be simplified. These are:
   * - If the object of type `T` is relatively small (less than 256 bytes) and
   *   if `T` satisfies `std::is_trivially_copyable`, then it is copied bit
   *   by bit into the output buffer.
   * - If no compression is requested, and if the object is a vector of objects
   *   whose type `T` satisfies `std::is_trivially_copyable`, then packing
   *   implies copying the length of the vector into the destination buffer
   *   followed by a bit-by-bit copy of the contents of the vector. A
   *   similar process is used for vectors of vectors of objects whose type
   *   `T` satisfies `std::is_trivially_copyable`.
   * - Finally, if the type `T` of the object to be packed is std::tuple<>
   *   (i.e., a tuple without any elements as indicated by the empty argument
   *   list) and if no compression is requested, then this
   *   type is considered an "empty" type and it is packed
   *   into a zero byte buffer. Using empty types is occasionally useful when
   *   sending messages to other processes if the important part about the
   *   message is that it is *sent*, not what it *contains* -- in other words,
   *   it puts the receiver on notice of something, without having to provide
   *   any details. In such cases, it is helpful if the message body can be
   *   empty -- that is, have length zero -- and using std::tuple<> facilitates
   *   this by providing a type which the present function packs into an
   *   empty output buffer, given that many deal.II functions send objects
   *   only after calling pack() to serialize them.
   *
   * In several of the special cases above, the `std::is_trivially_copyable`
   * property is important, see
   * https://en.cppreference.com/w/cpp/types/is_trivially_copyable .
   * For a type `T` to satisfy this property essentially means that an object
   * `t2` of this type can be initialized by copying another object `t1`
   * bit-by-bit into the memory space of `t2`. In particular, this is the case
   * for built-in types such as `int`, `double`, or `char`, as well as
   * structures and classes that only consist of such types and that have
   * neither user-defined constructors nor `virtual` functions. In practice,
   * and together with the fact that vectors and vector-of-vectors of these
   * types are also special-cased, this covers many of the most common kinds of
   * messages one sends around with MPI or one wants to serialize (the two
   * most common use cases for this function).
   */
  template <typename T>
  std::size_t
  pack(const T           &object,
       std::vector<char> &dest_buffer,
       const bool         allow_compression = true);

  /**
   * Creates and returns a buffer solely for the given object, using the
   * above mentioned pack function (including all of its special cases).
   *
   * If the library has been compiled with ZLIB enabled, then the output buffer
   * can be compressed. This can be triggered with the parameter
   * @p allow_compression, and is only of effect if ZLIB is enabled.
   */
  template <typename T>
  std::vector<char>
  pack(const T &object, const bool allow_compression = true);

  /**
   * Given a vector of characters, obtained through a call to the function
   * Utilities::pack, restore its content in an object of type `T`.
   *
   * This function uses boost::serialization utilities to unpack the object
   * from a vector of characters, and it is the inverse of the function
   * Utilities::pack(). It considers the same set of special cases as
   * documented with the pack() function.
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
   */
  template <typename T, int N>
  void
  unpack(const std::vector<char>::const_iterator &cbegin,
         const std::vector<char>::const_iterator &cend,
         T (&unpacked_object)[N],
         const bool allow_compression = true);

  /**
   * Check if the bit at position @p n in @p number is set.
   */
  bool
  get_bit(const unsigned char number, const unsigned int n);


  /**
   * Set the bit at position @p n in @p number to value @p x.
   */
  void
  set_bit(unsigned char &number, const unsigned int n, const bool x);


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
  dynamic_unique_cast(std::unique_ptr<From> &&p);

  /**
   * Return underlying value. Default: return input.
   */
  template <typename T>
  T &
  get_underlying_value(T &p);

  /**
   * Return underlying value. Specialization for std::shared_ptr<T>.
   */
  template <typename T>
  T &
  get_underlying_value(std::shared_ptr<T> &p);

  /**
   * Return underlying value. Specialization for const std::shared_ptr<T>.
   */
  template <typename T>
  T &
  get_underlying_value(const std::shared_ptr<T> &p);

  /**
   * Return underlying value. Specialization for std::unique_ptr<T>.
   */
  template <typename T>
  T &
  get_underlying_value(std::unique_ptr<T> &p);

  /**
   * Return underlying value. Specialization for const std::unique_ptr<T>.
   */
  template <typename T>
  T &
  get_underlying_value(const std::unique_ptr<T> &p);

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
    std::string
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
} // namespace Utilities


// --------------------- inline functions

namespace Utilities
{
  template <int N, typename T>
  inline constexpr T
  fixed_power(const T x)
  {
    Assert(((std::is_integral_v<T> == true) && (N >= 0)) ||
             (std::is_integral_v<T> == false),
           ExcMessage("If the type of the argument, T, is an integer type, "
                      "then the exponent N must be a non-negative integer "
                      "because the result would otherwise not be an integer."));

    if (N == 0)
      return T(1.);
    else if (N < 0)
      // For negative exponents, turn things into a positive exponent
      return T(1.) / fixed_power<-N>(x);
    else
      // If we get here, we have a positive exponent. Compute the result
      // by repeated squaring:
      return ((N % 2 == 1) ? x * fixed_power<N / 2>(x * x) :
                             fixed_power<N / 2>(x * x));
  }



  template <typename T, typename>
  constexpr DEAL_II_HOST_DEVICE T
  pow(const T base, const int iexp)
  {
#if defined(DEBUG) && !defined(DEAL_II_CXX14_CONSTEXPR_BUG)
    // Up to __builtin_expect this is the same code as in the 'Assert' macro.
    // The call to __builtin_expect turns out to be problematic.
#  if DEAL_II_KOKKOS_VERSION_GTE(3, 6, 0)
    KOKKOS_IF_ON_HOST(({
      if (!(iexp >= 0))
        ::dealii::deal_II_exceptions::internals::issue_error_noreturn(
          ::dealii::deal_II_exceptions::internals::ExceptionHandling::
            abort_or_throw_on_exception,
          __FILE__,
          __LINE__,
          __PRETTY_FUNCTION__,
          "iexp>=0",
          "ExcMessage(\"The exponent must not be negative!\")",
          ExcMessage("The exponent must not be negative!"));
    }))
#  else
#    ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    if (!(iexp >= 0))
      ::dealii::deal_II_exceptions::internals::issue_error_noreturn(
        ::dealii::deal_II_exceptions::internals::ExceptionHandling::
          abort_or_throw_on_exception,
        __FILE__,
        __LINE__,
        __PRETTY_FUNCTION__,
        "iexp>=0",
        "ExcMessage(\"The exponent must not be negative!\")",
        ExcMessage("The exponent must not be negative!"));
#    endif
#  endif
#endif
    // The "exponentiation by squaring" algorithm used below has to be expressed
    // in an iterative version since SYCL doesn't allow recursive functions used
    // in device code.

    if (iexp <= 0)
      return 1;

    int exp = iexp;
    T   x   = base;
    T   y   = 1;
    while (exp > 1)
      {
        if (exp % 2 == 1)
          y *= x;
        x *= x;
        exp /= 2;
      }
    return x * y;
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
                  DEAL_II_ASSERT_UNREACHABLE();
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

  namespace internal
  {
    /**
     * A structure that is used to identify whether a template argument is a
     * std::vector<T> or std::vector<std::vector<T>> where T is a type that
     * satisfies std::is_trivially_copyable_v<T> == true.
     */
    template <typename T>
    struct IsVectorOfTriviallyCopyable
    {
      static constexpr bool value = false;
    };



    template <typename T>
    struct IsVectorOfTriviallyCopyable<std::vector<T>>
    {
      static constexpr bool value =
        std::is_trivially_copyable_v<T> && !std::is_same_v<T, bool>;
    };



    template <typename T>
    struct IsVectorOfTriviallyCopyable<std::vector<std::vector<T>>>
    {
      static constexpr bool value =
        std::is_trivially_copyable_v<T> && !std::is_same_v<T, bool>;
    };



    /**
     * A function that is used to append the contents of a std::vector<T>
     * (where T is a type that satisfies std::is_trivially_copyable_v<T>
     * == true but not T==bool) bit for bit to a character array.
     *
     * If the type is not such a vector of T, then the function
     * throws an exception.
     */
    template <typename T>
    inline void
    append_vector_of_trivially_copyable_to_buffer(const T &,
                                                  std::vector<char> &)
    {
      // We shouldn't get here:
      DEAL_II_ASSERT_UNREACHABLE();
    }



    template <typename T,
              typename = std::enable_if_t<!std::is_same_v<T, bool> &&
                                          std::is_trivially_copyable_v<T>>>
    inline void
    append_vector_of_trivially_copyable_to_buffer(
      const std::vector<T> &object,
      std::vector<char>    &dest_buffer)
    {
      const typename std::vector<T>::size_type vector_size = object.size();

      // Reserve for the buffer so that it can store the size of 'object' as
      // well as all of its elements.
      dest_buffer.reserve(dest_buffer.size() + sizeof(vector_size) +
                          vector_size * sizeof(T));

      // Copy the size into the vector
      dest_buffer.insert(dest_buffer.end(),
                         reinterpret_cast<const char *>(&vector_size),
                         reinterpret_cast<const char *>(&vector_size + 1));

      // Insert the elements at the end of the vector:
      if (vector_size > 0)
        dest_buffer.insert(dest_buffer.end(),
                           reinterpret_cast<const char *>(object.data()),
                           reinterpret_cast<const char *>(object.data() +
                                                          vector_size));
    }



    template <typename T,
              typename = std::enable_if_t<!std::is_same_v<T, bool> &&
                                          std::is_trivially_copyable_v<T>>>
    inline void
    append_vector_of_trivially_copyable_to_buffer(
      const std::vector<std::vector<T>> &object,
      std::vector<char>                 &dest_buffer)
    {
      using size_type             = typename std::vector<T>::size_type;
      const size_type vector_size = object.size();

      typename std::vector<T>::size_type aggregated_size = 0;
      std::vector<size_type>             sizes;
      sizes.reserve(vector_size);
      for (const auto &a : object)
        {
          aggregated_size += a.size();
          sizes.push_back(a.size());
        }

      // Reserve for the buffer so that it can store the size of 'object' as
      // well as all of its elements.
      dest_buffer.reserve(dest_buffer.size() +
                          sizeof(vector_size) * (1 + vector_size) +
                          aggregated_size * sizeof(T));

      // Copy the size into the vector
      dest_buffer.insert(dest_buffer.end(),
                         reinterpret_cast<const char *>(&vector_size),
                         reinterpret_cast<const char *>(&vector_size + 1));

      // Copy the sizes of the individual chunks into the vector
      if (vector_size > 0)
        dest_buffer.insert(dest_buffer.end(),
                           reinterpret_cast<const char *>(sizes.data()),
                           reinterpret_cast<const char *>(sizes.data() +
                                                          vector_size));

      // Insert the elements at the end of the vector:
      for (const auto &a : object)
        dest_buffer.insert(dest_buffer.end(),
                           reinterpret_cast<const char *>(a.data()),
                           reinterpret_cast<const char *>(a.data() + a.size()));
    }



    template <typename T>
    inline void
    create_vector_of_trivially_copyable_from_buffer(
      const std::vector<char>::const_iterator &,
      const std::vector<char>::const_iterator &,
      T &)
    {
      // We shouldn't get here:
      DEAL_II_ASSERT_UNREACHABLE();
    }



    template <typename T,
              typename = std::enable_if_t<!std::is_same_v<T, bool> &&
                                          std::is_trivially_copyable_v<T>>>
    inline void
    create_vector_of_trivially_copyable_from_buffer(
      const std::vector<char>::const_iterator &cbegin,
      const std::vector<char>::const_iterator &cend,
      std::vector<T>                          &object)
    {
      // The size of the object vector can be found in cbegin of the buffer.
      // The data starts at cbegin + sizeof(vector_size).

      // Get the size of the vector
      typename std::vector<T>::size_type vector_size;
      std::memcpy(&vector_size, &*cbegin, sizeof(vector_size));

      Assert(static_cast<std::ptrdiff_t>(cend - cbegin) ==
               static_cast<std::ptrdiff_t>(sizeof(vector_size) +
                                           vector_size * sizeof(T)),
             ExcMessage("The given buffer has the wrong size."));
      (void)cend;

      // Resize the output array and copy the elements into it. We need
      // to make sure that we don't access the buffer via a
      // reinterpret_cast<T*> because the T objects in the buffer may
      // not be aligned properly. Rather, use memcpy, which doesn't care
      // about alignment and is correct in this situation because we
      // are dealing with trivially copyable objects.
      //
      // (Strictly speaking, this writes into the output object twice, once
      // for the resize() operation and once during memcpy. This could be
      // avoided by (i) checking whether the point is aligned, using for
      // example boost::alignment::is_aligned(), and (ii) if the data
      // is aligned, re-create the output vector using
      //   object.clear();
      //   object.insert (object.end(),
      //                  (T*)(&*cbegin + sizeof(vector_size)),
      //                  (T*)(&*cbegin + sizeof(vector_size)) + vector_size);
      // In practice, the difference is likely rather small, assuming the
      // compiler does not already optimize away the first initialization.
      object.resize(vector_size);
      if (vector_size > 0)
        std::memcpy(object.data(),
                    &*cbegin + sizeof(vector_size),
                    vector_size * sizeof(T));
    }



    template <typename T,
              typename = std::enable_if_t<!std::is_same_v<T, bool> &&
                                          std::is_trivially_copyable_v<T>>>
    inline void
    create_vector_of_trivially_copyable_from_buffer(
      const std::vector<char>::const_iterator &cbegin,
      const std::vector<char>::const_iterator &cend,
      std::vector<std::vector<T>>             &object)
    {
      // First get the size of the vector, and resize the output object
      using size_type = typename std::vector<T>::size_type;
      std::vector<char>::const_iterator iterator = cbegin;
      size_type                         vector_size;
      std::memcpy(&vector_size, &*iterator, sizeof(vector_size));
      object.clear();
      object.resize(vector_size);
      std::vector<size_type> sizes(vector_size);
      if (vector_size > 0)
        std::memcpy(sizes.data(),
                    &*iterator + sizeof(vector_size),
                    vector_size * sizeof(size_type));

      iterator += sizeof(vector_size) * (1 + vector_size);
      size_type aggregated_size = 0;
      for (const auto a : sizes)
        aggregated_size += a;

      Assert(static_cast<std::ptrdiff_t>(cend - iterator) ==
               static_cast<std::ptrdiff_t>(aggregated_size * sizeof(T)),
             ExcMessage("The given buffer has the wrong size."));
      (void)cend;

      // Then copy the elements. As for the previous function, use memcpy
      // rather than accessing the data via reinterpret_cast<T*> to
      // avoid alignment issues.
      for (unsigned int i = 0; i < vector_size; ++i)
        if (sizes[i] > 0)
          {
            object[i].resize(sizes[i]);
            std::memcpy(object[i].data(), &*iterator, sizes[i] * sizeof(T));
            iterator += sizes[i] * sizeof(T);
          }

      Assert(iterator == cend,
             ExcMessage("The given buffer has the wrong size."));
    }

  } // namespace internal



  template <typename T>
  std::size_t
  pack(const T           &object,
       std::vector<char> &dest_buffer,
       const bool         allow_compression)
  {
    std::size_t size = 0;


    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
    if constexpr (std::is_trivially_copyable<T>() && sizeof(T) < 256)
      {
        // Determine the size. There are places where we would like to use a
        // truly empty type, for which we use std::tuple<> (i.e., a tuple
        // of zero elements). For this class, the compiler reports a nonzero
        // sizeof(...) because that is the minimum possible for objects --
        // objects need to have distinct addresses, so they need to have a size
        // of at least one. But we can special case this situation.
        size = (std::is_same_v<T, std::tuple<>> ? 0 : sizeof(T));

        (void)allow_compression;
        const std::size_t previous_size = dest_buffer.size();
        dest_buffer.resize(previous_size + size);

        if (size > 0)
          std::memcpy(dest_buffer.data() + previous_size, &object, size);
      }
    // Next try if we have a vector of trivially copyable objects.
    // If that is the case, we can shortcut the whole BOOST serialization
    // machinery and just copy the content of the vector bit for bit
    // into the output buffer, assuming that we are not asked to compress
    // the data.
    else if (internal::IsVectorOfTriviallyCopyable<T>::value &&
             (allow_compression == false))
      {
        const std::size_t previous_size = dest_buffer.size();

        // When we have DEAL_II_HAVE_CXX17 set by default, we can just
        // inline the code of the following function here and make the 'if'
        // above a 'if constexpr'. Without the 'constexpr', we need to keep
        // the general template of the function that throws an exception.
        internal::append_vector_of_trivially_copyable_to_buffer(object,
                                                                dest_buffer);

        size = dest_buffer.size() - previous_size;
      }
    else
      {
        // use buffer as the target of a compressing
        // stream into which we serialize the current object
        const std::size_t previous_size = dest_buffer.size();
        {
          boost::iostreams::filtering_ostreambuf fosb;
#ifdef DEAL_II_WITH_ZLIB
          if (allow_compression)
            fosb.push(boost::iostreams::gzip_compressor());
#else
          (void)allow_compression;
#endif
          fosb.push(boost::iostreams::back_inserter(dest_buffer));

          boost::archive::binary_oarchive boa(fosb);
          boa << object;
          // the stream object has to be destroyed before the return statement
          // to ensure that all data has been written in the buffer
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
    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
    if constexpr (std::is_trivially_copyable<T>() && sizeof(T) < 256)
      {
        // Determine the size. There are places where we would like to use a
        // truly empty type, for which we use std::tuple<> (i.e., a tuple
        // of zero elements). For this class, the compiler reports a nonzero
        // sizeof(...) because that is the minimum possible for objects --
        // objects need to have distinct addresses, so they need to have a size
        // of at least one. But we can special case this situation.
        const std::size_t size =
          (std::is_same_v<T, std::tuple<>> ? 0 : sizeof(T));

        T object;

        (void)allow_compression;
        Assert(std::distance(cbegin, cend) == size, ExcInternalError());

        if (size > 0)
          std::memcpy(&object, &*cbegin, size);

        return object;
      }
    // Next try if we have a vector of trivially copyable objects.
    // If that is the case, we can shortcut the whole BOOST serialization
    // machinery and just copy the content of the buffer bit for bit
    // into an appropriately sized output vector, assuming that we
    // are not asked to compress the data.
    else if (internal::IsVectorOfTriviallyCopyable<T>::value &&
             (allow_compression == false))
      {
        // When we have DEAL_II_HAVE_CXX17 set by default, we can just
        // inline the code of the following function here and make the 'if'
        // above a 'if constexpr'. Without the 'constexpr', we need to keep
        // the general template of the function that throws an exception.
        T object;
        internal::create_vector_of_trivially_copyable_from_buffer(cbegin,
                                                                  cend,
                                                                  object);
        return object;
      }
    else
      {
        // decompress the buffer section into the object
        boost::iostreams::filtering_istreambuf fisb;
#ifdef DEAL_II_WITH_ZLIB
        if (allow_compression)
          fisb.push(boost::iostreams::gzip_decompressor());
#else
        (void)allow_compression;
#endif
        fisb.push(boost::iostreams::array_source(&*cbegin, cend - cbegin));

        boost::archive::binary_iarchive bia(fisb);

        T object;
        bia >> object;
        return object;
      }

    return T();
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
    if constexpr (std::is_trivially_copyable<T>() && sizeof(T) * N < 256)
      {
        Assert(std::distance(cbegin, cend) == sizeof(T) * N,
               ExcInternalError());
        std::memcpy(unpacked_object, &*cbegin, sizeof(T) * N);
      }
    else
      {
        // decompress the buffer section into the object
        boost::iostreams::filtering_istreambuf fisb;
#ifdef DEAL_II_WITH_ZLIB
        if (allow_compression)
          fisb.push(boost::iostreams::gzip_decompressor());
#else
        (void)allow_compression;
#endif
        fisb.push(boost::iostreams::array_source(&*cbegin, cend - cbegin));

        boost::archive::binary_iarchive bia(fisb);
        bia >> unpacked_object;
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



  inline bool
  get_bit(const unsigned char number, const unsigned int n)
  {
    AssertIndexRange(n, 8);

    // source:
    // https://stackoverflow.com/questions/47981/how-do-you-set-clear-and-toggle-a-single-bit
    // "Checking a bit"
    return ((number >> n) & 1U) != 0u;
  }



  inline void
  set_bit(unsigned char &number, const unsigned int n, const bool x)
  {
    AssertIndexRange(n, 8);

    // source:
    // https://stackoverflow.com/questions/47981/how-do-you-set-clear-and-toggle-a-single-bit
    // "Changing the nth bit to x"
    number ^= (-static_cast<unsigned char>(x) ^ number) & (1U << n);
  }



  template <typename To, typename From>
  inline std::unique_ptr<To>
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



  template <typename T>
  inline T &
  get_underlying_value(T &p)
  {
    return p;
  }



  template <typename T>
  inline T &
  get_underlying_value(std::shared_ptr<T> &p)
  {
    return *p;
  }



  template <typename T>
  inline T &
  get_underlying_value(const std::shared_ptr<T> &p)
  {
    return *p;
  }



  template <typename T>
  inline T &
  get_underlying_value(std::unique_ptr<T> &p)
  {
    return *p;
  }



  template <typename T>
  inline T &
  get_underlying_value(const std::unique_ptr<T> &p)
  {
    return *p;
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
    // automatically.
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
