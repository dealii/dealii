//----------------------------  grid_out.all_dimensions.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_out.all_dimensions.cc  ---------------------------


#include <grid/grid_out.h>


namespace GridOutFlags
{
  DX::DX (const bool write_faces,
	  const bool write_all_faces) :
    write_faces (write_faces),
    write_all_faces (write_all_faces)
  {};

  
  
  Ucd::Ucd (const bool write_preamble,
	    const bool write_faces) :
		  write_preamble (write_preamble),
		  write_faces (write_faces)
  {};

  
  
  Gnuplot::Gnuplot (const bool write_cell_numbers,
		    const unsigned int n_boundary_face_points) :
		  write_cell_numbers (write_cell_numbers),
		  n_boundary_face_points(n_boundary_face_points)
  {};

  
  
  EpsFlagsBase::EpsFlagsBase (const SizeType     size_type,
			      const unsigned int size,
			      const double       line_width,
			      const bool color_lines_on_user_flag,
			      const unsigned int n_boundary_face_points) :
		  size_type (size_type),
		  size (size),
		  line_width (line_width),
		  color_lines_on_user_flag(color_lines_on_user_flag),
		  n_boundary_face_points(n_boundary_face_points)
  {};
  

  
  Eps<1>::Eps (const SizeType     size_type,
	       const unsigned int size,
	       const double       line_width,
	       const bool color_lines_on_user_flag,
	       const unsigned int n_boundary_face_points)
		  :
		  EpsFlagsBase(size_type, size, line_width,
			       color_lines_on_user_flag,
			       n_boundary_face_points)
  {};

  
  
  Eps<2>::Eps (const SizeType     size_type,
	       const unsigned int size,
	       const double       line_width,
	       const bool color_lines_on_user_flag,
	       const unsigned int n_boundary_face_points)
		  :
		  EpsFlagsBase(size_type, size, line_width,
			       color_lines_on_user_flag,
			       n_boundary_face_points)
  {};


  
  Eps<3>::Eps (const SizeType     size_type,
	       const unsigned int size,
	       const double       line_width,
	       const bool color_lines_on_user_flag,
	       const unsigned int n_boundary_face_points,
	       const double        azimut_angle,
	       const double        turn_angle)
		  :
		  EpsFlagsBase(size_type, size, line_width,
			       color_lines_on_user_flag,
			       n_boundary_face_points),
		  azimut_angle (azimut_angle),
		  turn_angle (turn_angle)
  {};
};  // end namespace GridOutFlags



void GridOut::set_flags (const GridOutFlags::DX &flags) 
{
  dx_flags = flags;
};



void GridOut::set_flags (const GridOutFlags::Ucd &flags) 
{
  ucd_flags = flags;
};



void GridOut::set_flags (const GridOutFlags::Gnuplot &flags) 
{
  gnuplot_flags = flags;
};



void GridOut::set_flags (const GridOutFlags::Eps<1> &flags) 
{
  eps_flags_1 = flags;
};



void GridOut::set_flags (const GridOutFlags::Eps<2> &flags) 
{
  eps_flags_2 = flags;
};



void GridOut::set_flags (const GridOutFlags::Eps<3> &flags) 
{
  eps_flags_3 = flags;
};



std::string
GridOut::default_suffix (const OutputFormat output_format) 
{
  switch (output_format) 
    {
      case dx:
        return ".dx";
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
GridOut::parse_output_format (const std::string &format_name)
{
  if (format_name == "dx")
    return dx;

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



std::string GridOut::get_output_format_names () 
{
  return "dx|gnuplot|eps|ucd";
};



unsigned int
GridOut::memory_consumption () const
{
  return (sizeof(ucd_flags) +
	  sizeof(gnuplot_flags) +
	  sizeof(eps_flags_1) +
	  sizeof(eps_flags_2) +
	  sizeof(eps_flags_3));
};
