// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


#include <deal.II/base/data_out_base.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/insert_linebreaks.hpp>
#include <boost/archive/iterators/ostream_iterator.hpp>
#include <boost/archive/iterators/transform_width.hpp>

#include <string>
#include <vector>

#include "../tests.h"

#include "patches.h"

// Output data on repetitions of the unit hypercube.
//
// This test does not actually write any DX data, which is a binary format:
// instead, we remove date and times from the file and then convert to
// base64. This is necessary since numdiff cannot read binary files (but can
// read base64).

void
remove_datetime_and_convert_to_base64(const std::string &in, std::ostream &out)
{
  // clear the date string. This line has variable width so we have to use
  // substrings to fully excise it
  const std::string date_line_start = "attribute \"created\" string";
  const std::string date_line_last  = "\n";
  const std::size_t date_pos_front  = in.find(date_line_start);
  AssertThrow(date_pos_front != std::string::npos, ExcInternalError());
  const std::size_t date_pos_back = in.find(date_line_last, date_pos_front);
  AssertThrow(date_pos_back != std::string::npos, ExcInternalError());

  std::ostringstream output_no_date_stream;
  std::copy(in.begin(),
            in.begin() + date_pos_front,
            std::ostream_iterator<char>(output_no_date_stream));
  std::copy(in.begin() + date_pos_back + 1,
            in.end(),
            std::ostream_iterator<char>(output_no_date_stream));
  const std::string output_not_b64 = output_no_date_stream.str();

  using namespace boost::archive::iterators;
  // convert to base64 with lines of length 72
  using base64_text =
    insert_linebreaks<base64_from_binary<transform_width<const char *, 6, 8>>,
                      72>;

  std::copy(base64_text(output_not_b64.c_str()),
            base64_text(output_not_b64.c_str() + output_not_b64.size()),
            std::ostream_iterator<char>(out));
  out << std::endl;
}


template <int dim, int spacedim>
void
check(DataOutBase::DXFlags flags, std::ostream &out)
{
  const unsigned int np = 4;

  std::vector<DataOutBase::Patch<dim, spacedim>> patches(np);

  create_patches(patches);

  std::vector<std::string> names(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";

  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vectors;

  std::ostringstream binary_out;
  DataOutBase::write_dx(patches, names, vectors, flags, binary_out);
  const std::string output = binary_out.str();

  remove_datetime_and_convert_to_base64(output, out);
}

template <int dim>
void
check_cont(unsigned int         ncells,
           unsigned int         nsub,
           DataOutBase::DXFlags flags,
           std::ostream &       out)
{
  std::vector<DataOutBase::Patch<dim, dim>> patches;

  create_continuous_patches(patches, ncells, nsub);

  std::vector<std::string> names(1);
  names[0] = "CutOff";

  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vectors;

  std::ostringstream binary_out;
  DataOutBase::write_dx(patches, names, vectors, flags, binary_out);
  std::string output = binary_out.str();

  remove_datetime_and_convert_to_base64(output, out);
}


template <int dim, int spacedim>
void
check_all(std::ostream &log)
{
  std::ostream &out = log;

  char                 name[100];
  const char *         format = "%d%d%s.dx";
  DataOutBase::DXFlags flags(false, false, false, false);
  if (dim == 2 && spacedim == 2)
    {
      sprintf(name, format, dim, spacedim, "ffffcont");
      out << "==============================\n"
          << name << "\n==============================\n";
      check_cont<dim>(4, 4, flags, out);
    }
  if (true)
    {
      sprintf(name, format, dim, spacedim, "ffff");
      out << "==============================\n"
          << name << "\n==============================\n";
      check<dim, spacedim>(flags, out);
    }
  flags.int_binary = true;
  if (true)
    {
      sprintf(name, format, dim, spacedim, "tfff");
      out << "==============================\n"
          << name << "\n==============================\n";
      check<dim, spacedim>(flags, out);
    }
  flags.coordinates_binary = true;
  if (true)
    {
      sprintf(name, format, dim, spacedim, "ttff");
      out << "==============================\n"
          << name << "\n==============================\n";
      check<dim, spacedim>(flags, out);
    }
  flags.data_binary = true;
  if (true)
    {
      sprintf(name, format, dim, spacedim, "tttf");
      out << "==============================\n"
          << name << "\n==============================\n";
      check<dim, spacedim>(flags, out);
    }
  flags.data_double = true;
  if (true)
    {
      sprintf(name, format, dim, spacedim, "tttf");
      out << "==============================\n"
          << name << "\n==============================\n";
      check<dim, spacedim>(flags, out);
    }
}

int
main()
{
  std::ofstream logfile("output");
  check_all<1, 1>(logfile);
  check_all<1, 2>(logfile);
  check_all<2, 2>(logfile);
  check_all<2, 3>(logfile);
  check_all<3, 3>(logfile);
}
