//----------------------------  grid_out.all_dimensions.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_out.all_dimensions.cc  ---------------------------


#include <grid/grid_out.h>


GridOut::UcdFlags::UcdFlags (const bool write_preamble,
			     const bool write_faces) :
		write_preamble (write_preamble),
		write_faces (write_faces)
{};


GridOut::GnuplotFlags::GnuplotFlags (const bool write_cell_numbers,
				     const bool write_level) :
		write_cell_numbers (write_cell_numbers),
		write_level (write_level)
{};


GridOut::EpsFlagsBase::EpsFlagsBase (const SizeType     size_type,
				     const unsigned int size,
				     const double       line_width) :
		size_type (size_type),
		size (size),
		line_width (line_width)
{};


// egcs 1.1.2 does not understand this, so I made it inlined in the
// class declaration
//
// template <>
// GridOut::EpsFlags<3>::EpsFlags (const double        azimut_angle,
// 				const double        turn_angle) :
// 		azimut_angle (azimut_angle),
// 		turn_angle (turn_angle)
// {};


void GridOut::set_flags (const UcdFlags &flags) 
{
  ucd_flags = flags;
};


void GridOut::set_flags (const GnuplotFlags &flags) 
{
  gnuplot_flags = flags;
};


void GridOut::set_flags (const EpsFlags<1> &flags) 
{
  eps_flags_1 = flags;
};


void GridOut::set_flags (const EpsFlags<2> &flags) 
{
  eps_flags_2 = flags;
};


void GridOut::set_flags (const EpsFlags<3> &flags) 
{
  eps_flags_3 = flags;
};


string GridOut::default_suffix (const OutputFormat output_format) 
{
  switch (output_format) 
    {
      case gnuplot: 
	    return ".gnuplot";

      case ucd: 
	    return ".inp";

      case eps: 
	    return ".eps";
	    
      default: 
	    Assert (false, ExcNotImplemented()); 
	    return "";
    };
};


GridOut::OutputFormat
GridOut::parse_output_format (const string &format_name)
{
  if (format_name == "ucd")
    return ucd;

  if (format_name == "gnuplot")
    return gnuplot;

  if (format_name == "eps")
    return eps;
  
  AssertThrow (false, ExcInvalidState ());
				   // return something weird
  return OutputFormat(-1);
};


string GridOut::get_output_format_names () {
  return "gnuplot|eps";
};


