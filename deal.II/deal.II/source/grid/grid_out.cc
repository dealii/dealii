/* $Id$ */


#include <base/point.h>
#include <basic/grid_out.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>

#include <iomanip>
#include <algorithm>
#include <list>
#include <ctime>
#include <cmath>



template <int dim>
void GridOut::write_gnuplot (const Triangulation<dim> &tria,
			     ostream                  &out) 
{
  AssertThrow (out, ExcIO());

  typename Triangulation<dim>::active_cell_iterator        cell=tria.begin_active();
  const typename Triangulation<dim>::active_cell_iterator  endc=tria.end();
  for (; cell!=endc; ++cell)
    {
      out << "# cell " << cell << endl;
      
      switch (dim)
	{
	  case 1:
		out << cell->vertex(0) << ' ' << cell->level() << endl
		    << cell->vertex(1) << ' ' << cell->level() << endl
		    << endl;
		break;

	  case 2:
		out << cell->vertex(0) << ' ' << cell->level() << endl
		    << cell->vertex(1) << ' ' << cell->level() << endl
		    << cell->vertex(2) << ' ' << cell->level() << endl
		    << cell->vertex(3) << ' ' << cell->level() << endl
		    << cell->vertex(0) << ' ' << cell->level() << endl
		    << endl  // double new line for gnuplot 3d plots
		    << endl;
		break;

	  case 3:
						 // front face
		out << cell->vertex(0) << ' ' << cell->level() << endl
		    << cell->vertex(1) << ' ' << cell->level() << endl
		    << cell->vertex(2) << ' ' << cell->level() << endl
		    << cell->vertex(3) << ' ' << cell->level() << endl
		    << cell->vertex(0) << ' ' << cell->level() << endl
		    << endl;
						 // back face
		out << cell->vertex(4) << ' ' << cell->level() << endl
		    << cell->vertex(5) << ' ' << cell->level() << endl
		    << cell->vertex(6) << ' ' << cell->level() << endl
		    << cell->vertex(7) << ' ' << cell->level() << endl
		    << cell->vertex(4) << ' ' << cell->level() << endl
		    << endl;

						 // now for the four connecting lines
		out << cell->vertex(0) << ' ' << cell->level() << endl
		    << cell->vertex(4) << ' ' << cell->level() << endl
		    << endl;
		out << cell->vertex(1) << ' ' << cell->level() << endl
		    << cell->vertex(5) << ' ' << cell->level() << endl
		    << endl;
		out << cell->vertex(2) << ' ' << cell->level() << endl
		    << cell->vertex(6) << ' ' << cell->level() << endl
		    << endl;
		out << cell->vertex(3) << ' ' << cell->level() << endl
		    << cell->vertex(7) << ' ' << cell->level() << endl
		    << endl;
		break;
	};
    };
  
  AssertThrow (out, ExcIO());
};



template <int dim>
void GridOut::write_eps (const Triangulation<dim> &tria,
			 ostream                  &out) 
{
  Assert (false, ExcNotImplemented());
};



#if deal_II_dimension == 2

template <>
void GridOut::write_eps (const Triangulation<2> &tria,
			 ostream                &out) 
{
  typedef list<pair<Point<2>,Point<2> > > LineList;
  
  
  AssertThrow (out, ExcIO());

				   // make up a list of lines by which
				   // we will construct the triangulation
  LineList line_list;
  Triangulation<2>::active_line_iterator line   =tria.begin_active_line ();
  Triangulation<2>::active_line_iterator endline=tria.end_line ();

  for (; line!=endline; ++line)
    line_list.push_back (make_pair(line->vertex(0),
				   line->vertex(1)));

				   // find out minimum and maximum x and
				   // y coordinates to compute offsets
				   // and scaling factors
  double x_min = tria.begin_active_line()->vertex(0)(0);
  double x_max = x_min;
  double y_min = tria.begin_active_line()->vertex(0)(1);
  double y_max = y_min;

  for (LineList::const_iterator line=line_list.begin();
       line!=line_list.end(); ++line)
    {
      x_min = min (x_min, line->first(0));
      x_min = min (x_min, line->second(0));

      x_max = max (x_max, line->first(0));
      x_max = max (x_max, line->second(0));

      y_min = min (y_min, line->first(1));
      y_min = min (y_min, line->second(1));

      y_max = max (y_max, line->first(1));
      y_max = max (y_max, line->second(1));
    };

				   // scale in x-direction such that
				   // in the output 0 <= x <= 300.
				   // don't scale in y-direction to
				   // preserve the shape of the
				   // triangulation
  const double scale = 300. / (x_max - x_min);
  
  
				   // now write preamble
  if (true) 
    {
				       // block this to have local
				       // variables destroyed after
				       // use
      time_t  time1= time (0);
      tm     *time = localtime(&time1); 
      out << "%!PS-Adobe-2.0 EPSF-1.2" << endl
	  << "%%Title: deal.II Output" << endl
	  << "%%Creator: the deal.II library" << endl
	  << "%%Creation Date: " 
	  << time->tm_year+1900 << "/"
	  << time->tm_mon+1 << "/"
	  << time->tm_mday << " - "
	  << time->tm_hour << ":"
	  << setw(2) << time->tm_min << ":"
	  << setw(2) << time->tm_sec << endl
	  << "%%BoundingBox: "
					 // lower left corner
	  << "0 0 "
					 // upper right corner
	  << "300 " << static_cast<unsigned int>( (y_max-y_min) * scale )
	  << endl;

				       // define some abbreviations to keep
				       // the output small:
				       // m=move turtle to
				       // x=execute line stroke
      out << "/m {moveto} bind def" << endl
	  << "/x {lineto stroke} bind def" << endl;
      
      out << "%%EndProlog" << endl
	  << endl;

				       // set fine lines
      out << "0.5 setlinewidth" << endl;
    };

				   // now write the lines
  const Point<2> offset(x_min, y_min);
  
  for (LineList::const_iterator line=line_list.begin();
       line!=line_list.end(); ++line)
    out << (line->first  - offset) * scale << " m "
	<< (line->second - offset) * scale << " x\n";

  out << "showpage" << endl;
  
  AssertThrow (out, ExcIO());
};

#endif




template <int dim>
void GridOut::write (const Triangulation<dim> &tria,
		     ostream                  &out,
		     OutputFormat              output_format)
{
  switch (output_format)
    {
      case gnuplot:
	    write_gnuplot (tria, out);
	    return;

      case eps:
	    write_eps (tria, out);
	    return;
    };

  Assert (false, ExcInternalError());
};



string GridOut::default_suffix (const OutputFormat output_format) 
{
  switch (output_format) 
    {
      case gnuplot: 
	    return ".gnuplot";
	    
      case eps: 
	    return ".eps";
	    
      default: 
	    Assert (false, ExcNotImplemented()); 
	    return "";
    };
};
  



GridOut::OutputFormat
GridOut::parse_output_format (const string &format_name) {
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



// explicit instantiations
template void GridOut::write_gnuplot (const Triangulation<deal_II_dimension> &, ostream &);
//template void GridOut::write_eps (const Triangulation<deal_II_dimension> &, ostream &);
template void GridOut::write (const Triangulation<deal_II_dimension> &, ostream &, OutputFormat);
