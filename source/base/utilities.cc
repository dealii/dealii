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

#include <deal.II/base/config.h>

// It's necessary to include winsock2.h before thread_local_storage.h,
// because Intel implementation of TBB includes winsock.h,
// and we'll get a conflict between winsock.h and winsock2.h otherwise.
#ifdef DEAL_II_MSVC
#  include <winsock2.h>
#endif

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/utilities.h>

#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/iostreams/copy.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <algorithm>
#include <bitset>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#if defined(DEAL_II_HAVE_UNISTD_H) && defined(DEAL_II_HAVE_GETHOSTNAME)
#  include <unistd.h>
#endif

#ifndef DEAL_II_MSVC
#  include <cstdlib>
#endif


#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_WITH_MPI
#    include <deal.II/lac/trilinos_parallel_block_vector.h>
#    include <deal.II/lac/trilinos_vector.h>
#    include <deal.II/lac/vector_memory.h>

#    include <Epetra_MpiComm.h>
#    include <Teuchos_DefaultComm.hpp>
#  endif
#  include <Epetra_SerialComm.h>
#  include <Teuchos_RCP.hpp>
#endif

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  DeclException2(ExcInvalidNumber2StringConversersion,
                 unsigned int,
                 unsigned int,
                 << "When trying to convert " << arg1 << " to a string with "
                 << arg2 << " digits");
  DeclException1(ExcInvalidNumber, unsigned int, << "Invalid number " << arg1);
  DeclException1(ExcCantConvertString,
                 std::string,
                 << "Can't convert the string " << arg1
                 << " to the desired type");


  std::string
  dealii_version_string()
  {
    return DEAL_II_PACKAGE_NAME " version " DEAL_II_PACKAGE_VERSION;
  }


  namespace
  {
    template <int dim,
              typename Number,
              int effective_dim,
              typename LongDouble,
              typename Integer>
    std::vector<std::array<std::uint64_t, effective_dim>>
    inverse_Hilbert_space_filling_curve_effective(
      const std::vector<Point<dim, Number>> &points,
      const Point<dim, Number> &             bl,
      const std::array<LongDouble, dim> &    extents,
      const std::bitset<dim> &               valid_extents,
      const int                              min_bits,
      const Integer                          max_int)
    {
      std::vector<std::array<Integer, effective_dim>> int_points(points.size());

      for (unsigned int i = 0; i < points.size(); ++i)
        {
          // convert into integers:
          unsigned int eff_d = 0;
          for (unsigned int d = 0; d < dim; ++d)
            if (valid_extents[d])
              {
                Assert(extents[d] > 0, ExcInternalError());
                const LongDouble v = (static_cast<LongDouble>(points[i][d]) -
                                      static_cast<LongDouble>(bl[d])) /
                                     extents[d];
                Assert(v >= 0. && v <= 1., ExcInternalError());
                AssertIndexRange(eff_d, effective_dim);
                int_points[i][eff_d] =
                  static_cast<Integer>(v * static_cast<LongDouble>(max_int));
                ++eff_d;
              }
        }

      // note that we call this with "min_bits"
      return inverse_Hilbert_space_filling_curve<effective_dim>(int_points,
                                                                min_bits);
    }
  } // namespace

  template <int dim, typename Number>
  std::vector<std::array<std::uint64_t, dim>>
  inverse_Hilbert_space_filling_curve(
    const std::vector<Point<dim, Number>> &points,
    const int                              bits_per_dim)
  {
    using Integer = std::uint64_t;
    // take floating point number hopefully with mantissa >= 64bit
    using LongDouble = long double;

    // return if there is nothing to do
    if (points.size() == 0)
      return std::vector<std::array<std::uint64_t, dim>>();

    // get bounding box:
    Point<dim, Number> bl = points[0], tr = points[0];
    for (const auto &p : points)
      for (unsigned int d = 0; d < dim; ++d)
        {
          const double cid = p[d];
          bl[d]            = std::min(cid, bl[d]);
          tr[d]            = std::max(cid, tr[d]);
        }

    std::array<LongDouble, dim> extents;
    std::bitset<dim>            valid_extents;
    for (unsigned int i = 0; i < dim; ++i)
      {
        extents[i] =
          static_cast<LongDouble>(tr[i]) - static_cast<LongDouble>(bl[i]);
        valid_extents[i] = (extents[i] > 0.);
      }

    // make sure our conversion from fractional coordinates to
    // Integers work as expected, namely our cast (LongDouble)max_int
    const int min_bits =
      std::min(bits_per_dim,
               std::min(std::numeric_limits<Integer>::digits,
                        std::numeric_limits<LongDouble>::digits));

    // based on that get the maximum integer:
    const Integer max_int = (min_bits == std::numeric_limits<Integer>::digits ?
                               std::numeric_limits<Integer>::max() :
                               (Integer(1) << min_bits) - 1);

    const unsigned int effective_dim = valid_extents.count();
    if (effective_dim == dim)
      {
        return inverse_Hilbert_space_filling_curve_effective<dim,
                                                             Number,
                                                             dim,
                                                             LongDouble,
                                                             Integer>(
          points, bl, extents, valid_extents, min_bits, max_int);
      }

    // various degenerate cases
    std::array<std::uint64_t, dim> zero_ind;
    for (unsigned int d = 0; d < dim; ++d)
      zero_ind[d] = 0;

    std::vector<std::array<std::uint64_t, dim>> ind(points.size(), zero_ind);
    // manually check effective_dim == 1 and effective_dim == 2
    if (dim == 3 && effective_dim == 2)
      {
        const auto ind2 =
          inverse_Hilbert_space_filling_curve_effective<dim,
                                                        Number,
                                                        2,
                                                        LongDouble,
                                                        Integer>(
            points, bl, extents, valid_extents, min_bits, max_int);

        for (unsigned int i = 0; i < ind.size(); ++i)
          for (unsigned int d = 0; d < 2; ++d)
            ind[i][d + 1] = ind2[i][d];

        return ind;
      }
    else if (effective_dim == 1)
      {
        const auto ind1 =
          inverse_Hilbert_space_filling_curve_effective<dim,
                                                        Number,
                                                        1,
                                                        LongDouble,
                                                        Integer>(
            points, bl, extents, valid_extents, min_bits, max_int);

        for (unsigned int i = 0; i < ind.size(); ++i)
          ind[i][dim - 1] = ind1[i][0];

        return ind;
      }

    // we should get here only if effective_dim == 0
    Assert(effective_dim == 0, ExcInternalError());

    // if the bounding box is degenerate in all dimensions,
    // can't do much but exit gracefully by setting index according
    // to the index of each point so that there is no re-ordering
    for (unsigned int i = 0; i < points.size(); ++i)
      ind[i][dim - 1] = i;

    return ind;
  }



  template <int dim>
  std::vector<std::array<std::uint64_t, dim>>
  inverse_Hilbert_space_filling_curve(
    const std::vector<std::array<std::uint64_t, dim>> &points,
    const int                                          bits_per_dim)
  {
    using Integer = std::uint64_t;

    std::vector<std::array<Integer, dim>> int_points(points);

    std::vector<std::array<Integer, dim>> res(int_points.size());

    // follow
    // J. Skilling, Programming the Hilbert curve, AIP Conf. Proc. 707, 381
    // (2004); http://dx.doi.org/10.1063/1.1751381 also see
    // https://stackoverflow.com/questions/499166/mapping-n-dimensional-value-to-a-point-on-hilbert-curve
    // https://gitlab.com/octopus-code/octopus/blob/develop/src/grid/hilbert.c
    // https://github.com/trilinos/Trilinos/blob/master/packages/zoltan/src/hsfc/hsfc_hilbert.c
    // (Zoltan_HSFC_InvHilbertXd)
    // https://github.com/aditi137/Hilbert/blob/master/Hilbert/hilbert.cpp

    // now we can map to 1D coordinate stored in Transpose format
    // adopt AxestoTranspose function from the paper, that
    // transforms in-place between geometrical axes and Hilbert transpose.
    // Example:   b=5 bits for each of n=3 coordinates.
    //            15-bit Hilbert integer = A B C D E F G H I J K L M N O is
    //            stored as its Transpose
    //                   X[0] = A D G J M                X[2]|
    //                   X[1] = B E H K N    <------->       | /X[1]
    //                   X[2] = C F I L O               axes |/
    //                          high  low                    0------ X[0]

    // Depth of the Hilbert curve
    Assert(bits_per_dim <= std::numeric_limits<Integer>::digits,
           ExcMessage("This integer type can not hold " +
                      std::to_string(bits_per_dim) + " bits."));

    const Integer M = Integer(1) << (bits_per_dim - 1); // largest bit

    for (unsigned int index = 0; index < int_points.size(); ++index)
      {
        auto &X = int_points[index];
        auto &L = res[index];

        // Inverse undo
        for (Integer q = M; q > 1; q >>= 1)
          {
            const Integer p = q - 1;
            for (unsigned int i = 0; i < dim; i++)
              {
                // invert
                if (X[i] & q)
                  {
                    X[0] ^= p;
                  }
                // exchange
                else
                  {
                    const Integer t = (X[0] ^ X[i]) & p;
                    X[0] ^= t;
                    X[i] ^= t;
                  }
              }
          }

        // Gray encode (inverse of decode)
        for (unsigned int i = 1; i < dim; i++)
          X[i] ^= X[i - 1];

        Integer t = 0;
        for (Integer q = M; q > 1; q >>= 1)
          if (X[dim - 1] & q)
            t ^= q - 1;
        for (unsigned int i = 0; i < dim; i++)
          X[i] ^= t;

        // now we need to go from index stored in transpose format to
        // consecutive format, which is better suited for comparators.
        // we could interleave into some big unsigned int...
        // https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
        // https://stackoverflow.com/questions/4431522/given-2-16-bit-ints-can-i-interleave-those-bits-to-form-a-single-32-bit-int
        // ...but we would loose spatial resolution!

        // interleave using brute force, follow TransposetoLine from
        // https://github.com/aditi137/Hilbert/blob/master/Hilbert/hilbert.cpp
        {
          Integer      p = M;
          unsigned int j = 0;
          for (unsigned int i = 0; i < dim; ++i)
            {
              L[i] = 0;
              // go through bits using a mask q
              for (Integer q = M; q > 0; q >>= 1)
                {
                  if (X[j] & p)
                    L[i] |= q;
                  if (++j == dim)
                    {
                      j = 0;
                      p >>= 1;
                    }
                }
            }
        }

      } // end of the loop over points

    return res;
  }



  template <int dim>
  std::uint64_t
  pack_integers(const std::array<std::uint64_t, dim> &index,
                const int                             bits_per_dim)
  {
    using Integer = std::uint64_t;

    AssertIndexRange(bits_per_dim * dim, 65);
    Assert(bits_per_dim > 0, ExcMessage("bits_per_dim should be positive"));

    const Integer mask = (Integer(1) << bits_per_dim) - 1;

    Integer res = 0;
    for (unsigned int i = 0; i < dim; ++i)
      {
        // take bits_per_dim from each integer and shift them
        const Integer v = (mask & index[dim - 1 - i]) << (bits_per_dim * i);
        res |= v;
      }
    return res;
  }



  std::string
  compress(const std::string &input)
  {
#ifdef DEAL_II_WITH_ZLIB
    namespace bio = boost::iostreams;

    std::stringstream compressed;
    std::stringstream origin(input);

    bio::filtering_streambuf<bio::input> out;
    out.push(bio::gzip_compressor(
      bio::gzip_params(boost::iostreams::gzip::default_compression)));
    out.push(origin);
    bio::copy(out, compressed);

    return compressed.str();
#else
    return input;
#endif
  }



  std::string
  decompress(const std::string &compressed_input)
  {
#ifdef DEAL_II_WITH_ZLIB
    namespace bio = boost::iostreams;

    std::stringstream compressed(compressed_input);
    std::stringstream decompressed;

    bio::filtering_streambuf<bio::input> out;
    out.push(bio::gzip_decompressor());
    out.push(compressed);
    bio::copy(out, decompressed);

    return decompressed.str();
#else
    return compressed_input;
#endif
  }



  std::string
  encode_base64(const std::vector<unsigned char> &binary_input)
  {
    using namespace boost::archive::iterators;
    using It = base64_from_binary<
      transform_width<std::vector<unsigned char>::const_iterator, 6, 8>>;
    auto base64 = std::string(It(binary_input.begin()), It(binary_input.end()));
    // Add padding.
    return base64.append((3 - binary_input.size() % 3) % 3, '=');
  }



  std::vector<unsigned char>
  decode_base64(const std::string &base64_input)
  {
    using namespace boost::archive::iterators;
    using It =
      transform_width<binary_from_base64<std::string::const_iterator>, 8, 6>;
    auto binary = std::vector<unsigned char>(It(base64_input.begin()),
                                             It(base64_input.end()));
    // Remove padding.
    auto length = base64_input.size();
    if (binary.size() > 2 && base64_input[length - 1] == '=' &&
        base64_input[length - 2] == '=')
      {
        binary.erase(binary.end() - 2, binary.end());
      }
    else if (binary.size() > 1 && base64_input[length - 1] == '=')
      {
        binary.erase(binary.end() - 1, binary.end());
      }
    return binary;
  }



  std::string
  int_to_string(const unsigned int value, const unsigned int digits)
  {
    return to_string(value, digits);
  }



  template <typename number>
  std::string
  to_string(const number value, const unsigned int digits)
  {
    // For integer data types, use the standard std::to_string()
    // function. On the other hand, that function is defined in terms
    // of std::sprintf, which does not use the usual std::iostream
    // interface and tries to render floating point numbers in awkward
    // ways (see
    // https://en.cppreference.com/w/cpp/string/basic_string/to_string). So
    // resort to boost::lexical_cast for all other types (in
    // particular for floating point types.
    std::string lc_string = (std::is_integral<number>::value ?
                               std::to_string(value) :
                               boost::lexical_cast<std::string>(value));

    if ((digits != numbers::invalid_unsigned_int) &&
        (lc_string.size() < digits))
      {
        // We have to add the padding zeroes in front of the number
        const unsigned int padding_position = (lc_string[0] == '-') ? 1 : 0;

        const std::string padding(digits - lc_string.size(), '0');
        lc_string.insert(padding_position, padding);
      }

    return lc_string;
  }



  std::string
  replace_in_string(const std::string &input,
                    const std::string &from,
                    const std::string &to)
  {
    if (from.empty())
      return input;

    std::string            out = input;
    std::string::size_type pos = out.find(from);

    while (pos != std::string::npos)
      {
        out.replace(pos, from.size(), to);
        pos = out.find(from, pos + to.size());
      }
    return out;
  }

  std::string
  trim(const std::string &input)
  {
    std::string::size_type left  = 0;
    std::string::size_type right = input.size() > 0 ? input.size() - 1 : 0;

    for (; left < input.size(); ++left)
      {
        if (!std::isspace(input[left]))
          {
            break;
          }
      }

    for (; right >= left; --right)
      {
        if (!std::isspace(input[right]))
          {
            break;
          }
      }

    return std::string(input, left, right - left + 1);
  }



  std::string
  dim_string(const int dim, const int spacedim)
  {
    if (dim == spacedim)
      return int_to_string(dim);
    else
      return int_to_string(dim) + "," + int_to_string(spacedim);
  }


  unsigned int
  needed_digits(const unsigned int max_number)
  {
    if (max_number > 0)
      return static_cast<int>(
        std::ceil(std::log10(std::fabs(max_number + 0.1))));

    return 1;
  }



  template <typename Number>
  Number
  truncate_to_n_digits(const Number number, const unsigned int n_digits)
  {
    AssertThrow(n_digits >= 1, ExcMessage("invalid parameter."));

    if (!(std::fabs(number) > std::numeric_limits<Number>::min()))
      return number;

    const int order =
      static_cast<int>(std::floor(std::log10(std::fabs(number))));

    const int shift = -order + static_cast<int>(n_digits) - 1;

    Assert(shift <= static_cast<int>(std::floor(
                      std::log10(std::numeric_limits<Number>::max()))),
           ExcMessage(
             "Overflow. Use a smaller value for n_digits and/or make sure "
             "that the absolute value of 'number' does not become too small."));

    const Number factor = std::pow(10.0, static_cast<Number>(shift));

    const Number number_cutoff = std::trunc(number * factor) / factor;

    return number_cutoff;
  }


  int
  string_to_int(const std::string &s_)
  {
    // trim whitespace on either side of the text if necessary
    std::string s = s_;
    while ((s.size() > 0) && (s[0] == ' '))
      s.erase(s.begin());
    while ((s.size() > 0) && (s[s.size() - 1] == ' '))
      s.erase(s.end() - 1);

    // Now convert and see whether we succeed. Note that strtol only
    // touches errno if an error occurred, so if we want to check
    // whether an error happened, we need to make sure that errno==0
    // before calling strtol since otherwise it may be that the
    // conversion succeeds and that errno remains at the value it
    // was before, whatever that was.
    char *p;
    errno       = 0;
    const int i = std::strtol(s.c_str(), &p, 10);

    // We have an error if one of the following conditions is true:
    // - strtol sets errno != 0
    // - The original string was empty (we could have checked that
    //   earlier already)
    // - The string has non-zero length and strtol converted the
    //   first part to something useful, but stopped converting short
    //   of the terminating '\0' character. This happens, for example,
    //   if the given string is "1234 abc".
    AssertThrow(!((errno != 0) || (s.size() == 0) ||
                  ((s.size() > 0) && (*p != '\0'))),
                ExcMessage("Can't convert <" + s + "> to an integer."));

    return i;
  }



  std::vector<int>
  string_to_int(const std::vector<std::string> &s)
  {
    std::vector<int> tmp(s.size());
    for (unsigned int i = 0; i < s.size(); ++i)
      tmp[i] = string_to_int(s[i]);
    return tmp;
  }



  double
  string_to_double(const std::string &s_)
  {
    // trim whitespace on either side of the text if necessary
    std::string s = s_;
    while ((s.size() > 0) && (s[0] == ' '))
      s.erase(s.begin());
    while ((s.size() > 0) && (s[s.size() - 1] == ' '))
      s.erase(s.end() - 1);

    // Now convert and see whether we succeed. Note that strtol only
    // touches errno if an error occurred, so if we want to check
    // whether an error happened, we need to make sure that errno==0
    // before calling strtol since otherwise it may be that the
    // conversion succeeds and that errno remains at the value it
    // was before, whatever that was.
    char *p;
    errno          = 0;
    const double d = std::strtod(s.c_str(), &p);

    // We have an error if one of the following conditions is true:
    // - strtod sets errno != 0
    // - The original string was empty (we could have checked that
    //   earlier already)
    // - The string has non-zero length and strtod converted the
    //   first part to something useful, but stopped converting short
    //   of the terminating '\0' character. This happens, for example,
    //   if the given string is "1.234 abc".
    AssertThrow(!((errno != 0) || (s.size() == 0) ||
                  ((s.size() > 0) && (*p != '\0'))),
                ExcMessage("Can't convert <" + s + "> to a double."));

    return d;
  }



  std::vector<double>
  string_to_double(const std::vector<std::string> &s)
  {
    std::vector<double> tmp(s.size());
    for (unsigned int i = 0; i < s.size(); ++i)
      tmp[i] = string_to_double(s[i]);
    return tmp;
  }



  std::vector<std::string>
  split_string_list(const std::string &s, const std::string &delimiter)
  {
    // keep the currently remaining part of the input string in 'tmp' and
    // keep chopping elements of the list off the front
    std::string tmp = s;

    // as discussed in the documentation, eat whitespace from the end
    // of the string
    while (tmp.length() != 0 && tmp[tmp.length() - 1] == ' ')
      tmp.erase(tmp.length() - 1, 1);

    // split the input list until it is empty. since in every iteration
    // 'tmp' is what's left of the string after the next delimiter,
    // and since we've stripped trailing space already, 'tmp' will
    // be empty at one point if 's' ended in a delimiter, even if
    // there was space after the last delimiter. this matches what's
    // discussed in the documentation
    std::vector<std::string> split_list;
    while (tmp.length() != 0)
      {
        std::string name;
        name = tmp;

        if (name.find(delimiter) != std::string::npos)
          {
            name.erase(name.find(delimiter), std::string::npos);
            tmp.erase(0, tmp.find(delimiter) + delimiter.size());
          }
        else
          tmp = "";

        // strip spaces from this element's front and end
        while ((name.length() != 0) && (name[0] == ' '))
          name.erase(0, 1);
        while (name.length() != 0 && name[name.length() - 1] == ' ')
          name.erase(name.length() - 1, 1);

        split_list.push_back(name);
      }

    return split_list;
  }


  std::vector<std::string>
  split_string_list(const std::string &s, const char delimiter)
  {
    std::string d = ",";
    d[0]          = delimiter;
    return split_string_list(s, d);
  }


  std::vector<std::string>
  break_text_into_lines(const std::string &original_text,
                        const unsigned int width,
                        const char         delimiter)
  {
    std::string              text = original_text;
    std::vector<std::string> lines;

    // remove trailing spaces
    while ((text.length() != 0) && (text[text.length() - 1] == delimiter))
      text.erase(text.length() - 1, 1);

    // then split the text into lines
    while (text.length() != 0)
      {
        // in each iteration, first remove
        // leading spaces
        while ((text.length() != 0) && (text[0] == delimiter))
          text.erase(0, 1);

        std::size_t pos_newline = text.find_first_of('\n', 0);
        if (pos_newline != std::string::npos && pos_newline <= width)
          {
            std::string line(text, 0, pos_newline);
            while ((line.length() != 0) &&
                   (line[line.length() - 1] == delimiter))
              line.erase(line.length() - 1, 1);
            lines.push_back(line);
            text.erase(0, pos_newline + 1);
            continue;
          }

        // if we can fit everything into one
        // line, then do so. otherwise, we have
        // to keep breaking
        if (text.length() < width)
          {
            // remove trailing spaces
            while ((text.length() != 0) &&
                   (text[text.length() - 1] == delimiter))
              text.erase(text.length() - 1, 1);
            lines.push_back(text);
            text = "";
          }
        else
          {
            // starting at position width, find the
            // location of the previous space, so
            // that we can break around there
            int location = std::min<int>(width, text.length() - 1);
            for (; location > 0; --location)
              if (text[location] == delimiter)
                break;

            // if there are no spaces, then try if
            // there are spaces coming up
            if (location == 0)
              for (location = std::min<int>(width, text.length() - 1);
                   location < static_cast<int>(text.length());
                   ++location)
                if (text[location] == delimiter)
                  break;

            // now take the text up to the found
            // location and put it into a single
            // line, and remove it from 'text'
            std::string line(text, 0, location);
            while ((line.length() != 0) &&
                   (line[line.length() - 1] == delimiter))
              line.erase(line.length() - 1, 1);
            lines.push_back(line);
            text.erase(0, location);
          }
      }

    return lines;
  }



  bool
  match_at_string_start(const std::string &name, const std::string &pattern)
  {
    if (pattern.size() > name.size())
      return false;

    for (unsigned int i = 0; i < pattern.size(); ++i)
      if (pattern[i] != name[i])
        return false;

    return true;
  }



  std::pair<int, unsigned int>
  get_integer_at_position(const std::string &name, const unsigned int position)
  {
    Assert(position < name.size(), ExcInternalError());

    const std::string test_string(name.begin() + position, name.end());

    std::istringstream str(test_string);

    int i;
    if (str >> i)
      {
        // compute the number of
        // digits of i. assuming it
        // is less than 8 is likely
        // ok
        if (i < 10)
          return std::make_pair(i, 1U);
        else if (i < 100)
          return std::make_pair(i, 2U);
        else if (i < 1000)
          return std::make_pair(i, 3U);
        else if (i < 10000)
          return std::make_pair(i, 4U);
        else if (i < 100000)
          return std::make_pair(i, 5U);
        else if (i < 1000000)
          return std::make_pair(i, 6U);
        else if (i < 10000000)
          return std::make_pair(i, 7U);
        else
          {
            Assert(false, ExcNotImplemented());
            return std::make_pair(-1, numbers::invalid_unsigned_int);
          }
      }
    else
      return std::make_pair(-1, numbers::invalid_unsigned_int);
  }



  double
  generate_normal_random_number(const double a, const double sigma)
  {
    // if no noise: return now
    if (sigma == 0)
      return a;

    // we would want to use rand(), but that function is not reentrant
    // in a thread context. one could use rand_r, but this does not
    // produce reproducible results between threads either (though at
    // least it is reentrant). these two approaches being
    // non-workable, use a thread-local random number generator here.
    // we could use std::mt19937 but doing so results in compiler-dependent
    // output.
    static Threads::ThreadLocalStorage<boost::mt19937> random_number_generator;
    return boost::normal_distribution<>(a,
                                        sigma)(random_number_generator.get());
  }



  namespace System
  {
#if defined(__linux__)

    double
    get_cpu_load()
    {
      std::ifstream cpuinfo;
      cpuinfo.open("/proc/loadavg");

      AssertThrow(cpuinfo, ExcIO());

      double load;
      cpuinfo >> load;

      return load;
    }

#else

    double
    get_cpu_load()
    {
      return 0.;
    }

#endif

    const std::string
    get_current_vectorization_level()
    {
      switch (DEAL_II_VECTORIZATION_WIDTH_IN_BITS)
        {
          case 0:
            return "disabled";
          case 128:
#ifdef __ALTIVEC__
            return "AltiVec";
#else
            return "SSE2";
#endif
          case 256:
            return "AVX";
          case 512:
            return "AVX512";
          default:
            AssertThrow(false,
                        ExcInternalError(
                          "Invalid DEAL_II_VECTORIZATION_WIDTH_IN_BITS."));
            return "ERROR";
        }
    }


    void
    get_memory_stats(MemoryStats &stats)
    {
      stats.VmPeak = stats.VmSize = stats.VmHWM = stats.VmRSS = 0;

      // parsing /proc/self/stat would be a
      // lot easier, but it does not contain
      // VmHWM, so we use /status instead.
#if defined(__linux__)
      std::ifstream file("/proc/self/status");
      std::string   line;
      std::string   name;
      while (!file.eof())
        {
          file >> name;
          if (name == "VmPeak:")
            file >> stats.VmPeak;
          else if (name == "VmSize:")
            file >> stats.VmSize;
          else if (name == "VmHWM:")
            file >> stats.VmHWM;
          else if (name == "VmRSS:")
            {
              file >> stats.VmRSS;
              break; // this is always the last entry
            }

          getline(file, line);
        }
#endif
    }



    std::string
    get_hostname()
    {
#if defined(DEAL_II_HAVE_UNISTD_H) && defined(DEAL_II_HAVE_GETHOSTNAME)
      const unsigned int N = 1024;
      char               hostname[N];
      gethostname(&(hostname[0]), N - 1);
#else
      std::string hostname("unknown");
#endif
      return hostname;
    }



    std::string
    get_time()
    {
      std::time_t time1 = std::time(nullptr);
      std::tm *   time  = std::localtime(&time1);

      std::ostringstream o;
      o << time->tm_hour << ":" << (time->tm_min < 10 ? "0" : "")
        << time->tm_min << ":" << (time->tm_sec < 10 ? "0" : "")
        << time->tm_sec;

      return o.str();
    }



    std::string
    get_date()
    {
      std::time_t time1 = std::time(nullptr);
      std::tm *   time  = std::localtime(&time1);

      std::ostringstream o;
      o << time->tm_year + 1900 << "/" << time->tm_mon + 1 << "/"
        << time->tm_mday;

      return o.str();
    }



    void
    posix_memalign(void **memptr, std::size_t alignment, std::size_t size)
    {
#ifndef DEAL_II_MSVC
      const int ierr = ::posix_memalign(memptr, alignment, size);

      AssertThrow(ierr == 0, ExcOutOfMemory(size));
      AssertThrow(*memptr != nullptr, ExcOutOfMemory(size));
#else
      // Windows does not appear to have posix_memalign. just use the
      // regular malloc in that case
      *memptr = malloc(size);
      (void)alignment;
      AssertThrow(*memptr != 0, ExcOutOfMemory(size));
#endif
    }



    bool
    job_supports_mpi()
    {
      return Utilities::MPI::job_supports_mpi();
    }
  } // namespace System


#ifdef DEAL_II_WITH_TRILINOS

  namespace Trilinos
  {
    const Epetra_Comm &
    comm_world()
    {
#  ifdef DEAL_II_WITH_MPI
      static Teuchos::RCP<Epetra_MpiComm> communicator =
        Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD), true);
#  else
      static Teuchos::RCP<Epetra_SerialComm> communicator =
        Teuchos::rcp(new Epetra_SerialComm(), true);
#  endif

      return *communicator;
    }



    const Teuchos::RCP<const Teuchos::Comm<int>> &
    tpetra_comm_self()
    {
#  ifdef DEAL_II_WITH_MPI
      static auto communicator = Teuchos::RCP<const Teuchos::Comm<int>>(
        new Teuchos::MpiComm<int>(MPI_COMM_SELF));
#  else
      static auto communicator =
        Teuchos::RCP<const Teuchos::Comm<int>>(new Teuchos::Comm<int>());
#  endif

      return communicator;
    }



    const Epetra_Comm &
    comm_self()
    {
#  ifdef DEAL_II_WITH_MPI
      static Teuchos::RCP<Epetra_MpiComm> communicator =
        Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_SELF), true);
#  else
      static Teuchos::RCP<Epetra_SerialComm> communicator =
        Teuchos::rcp(new Epetra_SerialComm(), true);
#  endif

      return *communicator;
    }



    Epetra_Comm *
    duplicate_communicator(const Epetra_Comm &communicator)
    {
#  ifdef DEAL_II_WITH_MPI

      // see if the communicator is in fact a
      // parallel MPI communicator; if so,
      // return a duplicate of it
      const Epetra_MpiComm *mpi_comm =
        dynamic_cast<const Epetra_MpiComm *>(&communicator);
      if (mpi_comm != nullptr)
        return new Epetra_MpiComm(
          Utilities::MPI::duplicate_communicator(mpi_comm->GetMpiComm()));
#  endif

      // if we don't support MPI, or if the
      // communicator in question was in fact
      // not an MPI communicator, return a
      // copy of the same object again
      Assert(dynamic_cast<const Epetra_SerialComm *>(&communicator) != nullptr,
             ExcInternalError());
      return new Epetra_SerialComm(
        dynamic_cast<const Epetra_SerialComm &>(communicator));
    }



    void
    destroy_communicator(Epetra_Comm &communicator)
    {
      // save the communicator, reset the map, and delete the communicator if
      // this whole thing was created as an MPI communicator
#  ifdef DEAL_II_WITH_MPI
      Epetra_MpiComm *mpi_comm = dynamic_cast<Epetra_MpiComm *>(&communicator);
      if (mpi_comm != nullptr)
        {
          MPI_Comm comm  = mpi_comm->GetMpiComm();
          *mpi_comm      = Epetra_MpiComm(MPI_COMM_SELF);
          const int ierr = MPI_Comm_free(&comm);
          AssertThrowMPI(ierr);
        }
#  endif
    }



    unsigned int
    get_n_mpi_processes(const Epetra_Comm &mpi_communicator)
    {
      return mpi_communicator.NumProc();
    }


    unsigned int
    get_this_mpi_process(const Epetra_Comm &mpi_communicator)
    {
      return static_cast<unsigned int>(mpi_communicator.MyPID());
    }



    Epetra_Map
    duplicate_map(const Epetra_BlockMap &map, const Epetra_Comm &comm)
    {
      if (map.LinearMap() == true)
        {
          // each processor stores a
          // contiguous range of
          // elements in the
          // following constructor
          // call
          return Epetra_Map(map.NumGlobalElements(),
                            map.NumMyElements(),
                            map.IndexBase(),
                            comm);
        }
      else
        {
          // the range is not
          // contiguous
          return Epetra_Map(map.NumGlobalElements(),
                            map.NumMyElements(),
                            map.MyGlobalElements(),
                            0,
                            comm);
        }
    }
  } // namespace Trilinos

#endif

  template std::string
  to_string<int>(int, unsigned int);
  template std::string
  to_string<long int>(long int, unsigned int);
  template std::string
  to_string<long long int>(long long int, unsigned int);
  template std::string
  to_string<unsigned int>(unsigned int, unsigned int);
  template std::string
  to_string<unsigned long int>(unsigned long int, unsigned int);
  template std::string
  to_string<unsigned long long int>(unsigned long long int, unsigned int);
  template std::string
  to_string<float>(float, unsigned int);
  template std::string
  to_string<double>(double, unsigned int);
  template std::string
  to_string<long double>(long double, unsigned int);

  template double
  truncate_to_n_digits(const double, const unsigned int);
  template float
  truncate_to_n_digits(const float, const unsigned int);

  template std::vector<std::array<std::uint64_t, 1>>
  inverse_Hilbert_space_filling_curve<1, double>(
    const std::vector<Point<1, double>> &,
    const int);
  template std::vector<std::array<std::uint64_t, 1>>
  inverse_Hilbert_space_filling_curve<1>(
    const std::vector<std::array<std::uint64_t, 1>> &,
    const int);
  template std::vector<std::array<std::uint64_t, 2>>
  inverse_Hilbert_space_filling_curve<2, double>(
    const std::vector<Point<2, double>> &,
    const int);
  template std::vector<std::array<std::uint64_t, 2>>
  inverse_Hilbert_space_filling_curve<2>(
    const std::vector<std::array<std::uint64_t, 2>> &,
    const int);
  template std::vector<std::array<std::uint64_t, 3>>
  inverse_Hilbert_space_filling_curve<3, double>(
    const std::vector<Point<3, double>> &,
    const int);
  template std::vector<std::array<std::uint64_t, 3>>
  inverse_Hilbert_space_filling_curve<3>(
    const std::vector<std::array<std::uint64_t, 3>> &,
    const int);

  template std::uint64_t
  pack_integers<1>(const std::array<std::uint64_t, 1> &, const int);
  template std::uint64_t
  pack_integers<2>(const std::array<std::uint64_t, 2> &, const int);
  template std::uint64_t
  pack_integers<3>(const std::array<std::uint64_t, 3> &, const int);


} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE
