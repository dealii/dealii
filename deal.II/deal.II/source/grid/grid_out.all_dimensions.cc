//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <grid/grid_out.h>


namespace GridOutFlags
{
  DX::DX (const bool write_cells,
	  const bool write_faces,
	  const bool write_diameter,
	  const bool write_measure,
	  const bool write_all_faces) :
    write_cells (write_cells),
    write_faces (write_faces),
    write_diameter (write_diameter),
    write_measure (write_measure),
    write_all_faces (write_all_faces)
  {}

  
  
  Ucd::Ucd (const bool write_preamble,
	    const bool write_faces) :
		  write_preamble (write_preamble),
		  write_faces (write_faces)
  {}

  
  
  Gnuplot::Gnuplot (const bool write_cell_numbers,
		    const unsigned int n_boundary_face_points) :
		  write_cell_numbers (write_cell_numbers),
		  n_boundary_face_points(n_boundary_face_points)
  {}

  
  
  EpsFlagsBase::EpsFlagsBase (const SizeType     size_type,
			      const unsigned int size,
			      const double       line_width,
			      const bool color_lines_on_user_flag,
                            const unsigned int n_boundary_face_points,
                            const bool color_lines_level) :
		  size_type (size_type),
		  size (size),
		  line_width (line_width),
		  color_lines_on_user_flag(color_lines_on_user_flag),
                n_boundary_face_points(n_boundary_face_points),
                color_lines_level(color_lines_level)
  {}
  

  
  Eps<1>::Eps (const SizeType     size_type,
	       const unsigned int size,
	       const double       line_width,
	       const bool color_lines_on_user_flag,
	       const unsigned int n_boundary_face_points)
		  :
		  EpsFlagsBase(size_type, size, line_width,
			       color_lines_on_user_flag,
			       n_boundary_face_points)
  {}

  
  
  Eps<2>::Eps (const SizeType     size_type,
	       const unsigned int size,
	       const double       line_width,
	       const bool color_lines_on_user_flag,
	       const unsigned int n_boundary_face_points,
	       const bool         write_cell_numbers,
	       const bool         write_cell_number_level,
             const bool         write_vertex_numbers,
             const bool         color_lines_level
      )
		  :
		  EpsFlagsBase(size_type, size, line_width,
			       color_lines_on_user_flag,
                             n_boundary_face_points,
                             color_lines_level),
		  write_cell_numbers (write_cell_numbers),
		  write_cell_number_level (write_cell_number_level),
		  write_vertex_numbers (write_vertex_numbers)
  {}


  
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
  {}


  XFig::XFig ()
		  :
    draw_boundary(true),
    level_color(false),
    level_depth(true),
    n_boundary_face_points(0),
    scaling(1.,1.),
    fill_style (20),
    line_style(0),
    line_thickness(1),
    boundary_style(0),
    boundary_thickness(3)
  {}
}  // end namespace GridOutFlags



void GridOut::set_flags (const GridOutFlags::DX &flags) 
{
  dx_flags = flags;
}



void GridOut::set_flags (const GridOutFlags::Ucd &flags) 
{
  ucd_flags = flags;
}



void GridOut::set_flags (const GridOutFlags::Gnuplot &flags) 
{
  gnuplot_flags = flags;
}



void GridOut::set_flags (const GridOutFlags::Eps<1> &flags) 
{
  eps_flags_1 = flags;
}



void GridOut::set_flags (const GridOutFlags::Eps<2> &flags) 
{
  eps_flags_2 = flags;
}



void GridOut::set_flags (const GridOutFlags::Eps<3> &flags) 
{
  eps_flags_3 = flags;
}



void GridOut::set_flags (const GridOutFlags::XFig &flags) 
{
  xfig_flags = flags;
}



std::string
GridOut::default_suffix (const OutputFormat output_format) 
{
  switch (output_format) 
    {
      case none:
	return "";
      case dx:
        return ".dx";
      case gnuplot: 
	return ".gnuplot";
      case ucd: 
	return ".inp";
      case eps: 
	return ".eps";
      case xfig:
	return ".fig";
      default: 
	    Assert (false, ExcNotImplemented()); 
	    return "";
    };
}



GridOut::OutputFormat
GridOut::parse_output_format (const std::string &format_name)
{
  if (format_name == "none" || format_name == "false")
    return none;

  if (format_name == "dx")
    return dx;

  if (format_name == "ucd")
    return ucd;

  if (format_name == "gnuplot")
    return gnuplot;

  if (format_name == "eps")
    return eps;
  
  if (format_name == "xfig")
    return xfig;
  
  AssertThrow (false, ExcInvalidState ());
				   // return something weird
  return OutputFormat(-1);
}



std::string GridOut::get_output_format_names () 
{
  return "none|dx|gnuplot|eps|ucd|xfig";
}



unsigned int
GridOut::memory_consumption () const
{
  return (sizeof(ucd_flags) +
	  sizeof(gnuplot_flags) +
	  sizeof(eps_flags_1) +
	  sizeof(eps_flags_2) +
	  sizeof(eps_flags_3) +
	  sizeof(xfig_flags));
}
