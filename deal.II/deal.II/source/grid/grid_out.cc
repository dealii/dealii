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
  typedef list<pair<Point<2>,Point<2> > > LineList;
  
  
  AssertThrow (out, ExcIO());

				   // make up a list of lines by which
				   // we will construct the triangulation
				   //
				   // this part unfortunately is a bit
				   // dimension dependent, so we have to
				   // treat every dimension different.
				   // however, by directly producing
				   // the lines to be printed, i.e. their
				   // 2d images, we can later do the
				   // actual output dimension independent
				   // again
  LineList line_list;

  switch (dim)
    {
      case 2:
      {
	Triangulation<dim>::active_line_iterator line   =tria.begin_active_line ();
	Triangulation<dim>::active_line_iterator endline=tria.end_line ();
	
	for (; line!=endline; ++line)
					   // one would expect
					   // make_pair(line->vertex(0),
					   //           line->vertex(1))
					   // here, but that is not
					   // dimension independent, since
					   // vertex(i) is Point<dim>,
					   // but we want a Point<2>.
					   // in fact, whenever we're here,
					   // the vertex is a Point<dim>,
					   // but the compiler does not
					   // know this. hopefully, the
					   // compiler will optimize away
					   // this little kludge
	  line_list.push_back (make_pair(Point<2>(line->vertex(0)(0),
						  line->vertex(0)(1)),
					 Point<2>(line->vertex(1)(0),
						  line->vertex(1)(1))));
	break;
      };
       
      case 3:
      {
	Triangulation<dim>::active_line_iterator line   =tria.begin_active_line ();
	Triangulation<dim>::active_line_iterator endline=tria.end_line ();
	
					 // loop over all lines and compute their
					 // projection on the plane perpendicular
					 // to the direction of sight

					 // direction of view equals the unit 
					 // vector of the position of the
					 // spectator to the origin.
					 //
					 // we chose here the viewpoint as in
					 // gnuplot
	const double z_angle    = 60;
	const double turn_angle = 30;
	const double pi = 3.1415926536;
	const Point<dim> view_direction(-sin(z_angle * 2.*pi / 360.) * sin(turn_angle * 2.*pi / 360.),
					+sin(z_angle * 2.*pi / 360.) * cos(turn_angle * 2.*pi / 360.),
					-cos(z_angle * 2.*pi / 360.));
	
					 // decide about the two unit vectors
					 // in this plane. we chose the first one
					 // to be the projection of the z-axis
					 // to this plane
	const Point<dim> vector1
	  = Point<dim>(0,0,1) - ((Point<dim>(0,0,1) * view_direction) * view_direction);
	const Point<dim> unit_vector1 = vector1 / sqrt(vector1.square());
	
					 // now the third vector is fixed. we
					 // chose the projection of a more or
					 // less arbitrary vector to the plane
					 // perpendicular to the first one
	const Point<dim> vector2
	  = (Point<dim>(1,0,0)
	     - ((Point<dim>(1,0,0) * view_direction) * view_direction)
	     - ((Point<dim>(1,0,0) * unit_vector1)   * unit_vector1));
	const Point<dim> unit_vector2 = vector2 / sqrt(vector2.square());
	
	for (; line!=endline; ++line) 
	  line_list.push_back (make_pair(Point<2>(line->vertex(0) * unit_vector2,
						  line->vertex(0) * unit_vector1),
					 Point<2>(line->vertex(1) * unit_vector2,
						  line->vertex(1) * unit_vector1)));

	break;
      };

      default:
	    Assert (false, ExcNotImplemented());
    };
  
  

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




// explicit instantiations
template void GridOut::write_gnuplot (const Triangulation<deal_II_dimension> &, ostream &);
template void GridOut::write_eps (const Triangulation<deal_II_dimension> &, ostream &);
template void GridOut::write (const Triangulation<deal_II_dimension> &, ostream &, OutputFormat);
