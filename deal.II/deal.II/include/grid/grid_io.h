/*----------------------------   grid_io.h     ---------------------------*/
/*      $Id$                 */
#ifndef __grid_io_H
#define __grid_io_H
/*----------------------------   grid_io.h     ---------------------------*/

#include <basic/forward-declarations.h>
#include <base/exceptions.h>
#include <string>



/**
 * This class provides a means to output a triangulation to a file in different
 * formats. Presently provided are functions to write it in GNUPLOT and
 * encapsulated postscript format. There are also function to dispatch to
 * the different output functions based on a parameter given, e.g., through
 * a parameter file, thus making user programs invariant of the number and
 * names of the file formats implemented in this class.
 *
 *
 * \subsection{GNUPLOT format}
 * In GNUPLOT format, each cell is written as a sequence of lines. There is not
 * much more to be said here since the actual representation as well as the
 * viewing angle (for 3d plots) is done by GNUPLOT itself.
 *
 * One additional feature, however, is worth being mentioned: after each point
 * denoting one of the end points of a line, the level of the respective cell
 * is written. Therefore, if you let GNUPLOT draw a 2d grid as a 3d plot, you
 * will see more refined cells being raised against cells with less refinement.
 * Also, if you draw a cut through a 3d grid, you can extrude the refinement
 * level in the direction orthogonal to the cut plane.
 *
 *
 * \subsection{Encapsulated postscript format}
 * In this format, each line of the triangulation is written separately. We
 * scale the picture such that x-values range between zero and 300 and the y-range
 * is scaled with the same factor. The bounding box is close to the triangulation
 * on all four sides, without an extra frame. The line width is chosen to be 0.5,
 * which is a thin line relative to the extension of the picture of 300.
 * 
 *
 * @author Wolfgang Bangerth, 1999; postscript format based on an implementation by Stefan Nauber, 1999
 */
class GridOut 
{
  public:
				     /**
				      * Declaration of a name for each of the
				      * different output formats.
				      */
    enum OutputFormat { gnuplot, eps };

				     /**
				      * Write the triangulation in the
				      * gnuplot format. See the general
				      * documentation for a description
				      * of what happens here.
				      */
    template <int dim>
    static void write_gnuplot (const Triangulation<dim> &tria,
			       ostream                  &out);

				     /**
				      * Write the triangulation in the
				      * encapsulated postscript format. See the
				      * general documentation for a description
				      * of what happens here.
				      */
    template <int dim>
    static void write_eps (const Triangulation<dim> &tria,
			   ostream                  &out);

				     /**
				      * Write data and grid to #out# according
				      * to the given data format. This function
				      * simply calls the appropriate
				      * #write_*# function.
				      */
    template <int dim>
    static void write (const Triangulation<dim> &tria,
		       ostream                  &out,
		       const OutputFormat        output_format);
    
				     /**
				      * Provide a function which tells us which
				      * suffix with a given output format
				      * usually has. At present the following
				      * formats are defined:
				      * \begin{itemize}
				      * \item #gnuplot#: #.gnuplot#
				      * \item #eps#: #.eps#.
				      * \end{itemize}
				      *
				      * Since this function does not need data
				      * from this object, it is static and can
				      * thus be called without creating an
				      * object of this class.
				      */
    static string default_suffix (const OutputFormat output_format);

				     /**
				      * Return the #OutputFormat# value
				      * corresponding to the given string. If
				      * the string does not match any known
				      * format, an exception is thrown.
				      *
				      * Since this function does not need data
				      * from this object, it is static and can
				      * thus be called without creating an
				      * object of this class. Its main purpose
				      * is to allow a program to use any
				      * implemented output format without the
				      * need to extend the program's parser
				      * each time a new format is implemented.
				      *
				      * To get a list of presently available
				      * format names, e.g. to give it to the
				      * #ParameterHandler# class, use the
				      * function #get_output_format_names ()#.
				      */
    static OutputFormat parse_output_format (const string &format_name);

				     /**
				      * Return a list of implemented output
				      * formats. The different names are
				      * separated by vertical bar signs (#`|'#)
				      * as used by the #ParameterHandler#
				      * classes.
				      */
    static string get_output_format_names ();


				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidState);
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);
};



/*----------------------------   grid_io.h     ---------------------------*/
/* end of #ifndef __grid_io_H */
#endif
/*----------------------------   grid_io.h     ---------------------------*/
