/*----------------------------   data_out_base.h     ---------------------------*/
/*      $Id$                 */
#ifndef __data_out_base_H
#define __data_out_base_H
/*----------------------------   data_out_base.h     ---------------------------*/


#include <base/point.h>
#include <lac/fullmatrix.h>
#include <grid/geometry_info.h>
#include <iostream>
#include <vector>
#include <string>






/**
 * This is a base class for output of data on meshes of very general
 * form. It basically only provides a set of functions for several output
 * formats which take a list of patches and write them to an output
 * stream.
 *
 * By offering this interface to the different output formats, it is simple
 * to extend this class to new formats without depending on such things
 * as actual triangulations and handling of data vectors. These things shall
 * be provided by derived class which have a user callable interface then.
 *
 *
 * \subsection{Interface}
 * This class has an interface that is not usually called by a user directly;
 * also, it consists of #static# functions only. Usually, derived classes will
 * inherit this class #protected# to hide this interface to the users of thes
 * classes.
 *
 * The interface of this class basically consists of the declaration of a data
 * type describing a patch and a bunch of functions taking a list of patches
 * and writing them in one format or other to the stream. It is in the
 * responsibility of the derived classes to provide this list of patches.
 * In addition to the list of patches, a name for each data set may be given.
 *
 *
 * \subsection{Patches}
 * Grids can be thought of as a collection of cells; if you want to write out
 * data on such a grid, you can do so by writing them one cell at a time.
 * The functions in this class therefore take a list of objects describing the
 * data on one cell each. This data for each cell usually consists of a list
 * of vertices for this cell, and a list of data values (for example solution
 * data, error information, etc) at each of these vertices.
 *
 * In some cases, this interface to a cell is too restricted, however. For
 * example, you may have higher order elements and printing the values at
 * the vertices only is not enough. For this reason, we not only provide
 * writing the data on the vertices only, but the data is organizes as a
 * tensor product grid on each cell. The parameter #n_subdivision#, which is
 * given for each patch separately, denotes how often the cell is to be
 * divided for output; for example, #n_subdivisions==1# yields no subdivision
 * of the cell, #n_subdivisions==2# will produce a grid of 3 times 3 points
 * in two spatial dimensions and 3 times 3 times 3 points in three dimensions,
 * #n_subdivisions==2# will yield 4 times 4 (times 4) points, etc. The actual
 * location of these points on the patch will be computed by a multilinear
 * transformation from the vertices given for this patch.
 *
 * Given these comments, the actual data to be printed on this patch of
 * points consists of several data sets each of which has a value at each
 * of the patch points. For example with #n_subdivisions==2# in two space
 * dimensions, each data set has to provide nine values, and since the
 * patch is to be printed as a tensor product (or its transformation to the
 * real space cell), its values are to be ordered like
 * #(x0,y0) (x0,y1) (x0,y2) (x1,y0) (x1,y1) (x1,y2) (x2,y0) (x2,y1) (x2,y2)#,
 * i.e. the z-coordinate runs fastest, then the y-coordinate, then x (if there
 * are that many space directions).
 *
 * The #Patch# class takes a template parameter denoting the space dimension
 * in which this patch operates.
 *
 *
 * \section{Supported output formats}
 *
 * \subsection{UCD format}
 *
 * The UCD format is described in the AVS developer's guide. Due to
 * limitations in the present format, only node based data can be output,
 * which in one reason why we invented the patch concept. In order to
 * write higher order elements, you may split them up into several subdivisions
 * of each cell. These subcells will then, however, also appear as different
 * cells by programs which understand the UCD format.
 * 
 * No use is made of the possibility to give model data since these
 * are not supported by all UCD aware programs. You may give cell data
 * in derived classes by setting all values of a given data set on a
 * patch to the same value.
 *
 *
 * \subsection{GNUPLOT format}
 *
 * The GNUPLOT format is not able to handle data on unstructured grids
 * directly. Directly would mean that you only give the vertices and
 * the solution values thereon and the program constructs its own grid
 * to represent the data. This is only possible for a structured tensor
 * product grid in two dimensions. However, it is possible to give several
 * such patches within one file, which is exactly what the respective
 * function of this class does: writing each cell's data as a patch of
 * data, at least if the patches as passed from derived classes
 * represent cells. Note that the functions on patches need not be
 * continuous at interfaces between patches, so this method also works
 * for discontinuous elements. Note also, that GNUPLOT can do hidden
 * line removal for patched data.
 *
 * While this discussion applies to two spatial dimensions, it is more
 * complicated in 3d. The reason is that we could still use patches, but
 * it is difficult when trying to visualize them, since if we use a cut
 * through the data (by, for example, using x- and z-coordinates, a fixed
 * y-value and plot function values in z-direction, then the patched data
 * is not a patch in the sense GNUPLOT wants it any more. Therefore, we use
 * another approach, namely writing the data on the 3d grid as a sequence
 * of lines, i.e. two points each associated with one or more data sets.
 * There are therefore 12 lines for each subcells of a patch.
 *
 * Given the lines as described above, a cut through this data in Gnuplot
 * can then be achieved like this:
 * \begin{verbatim}
 *   set data style lines
 *   splot [:][:][0:] "T" using 1:2:($3==.5 ? $4 : -1)
 * \end{verbatim}
 * This command plots data in x- and y-direction unbounded, but in z-direction
 * only those data points which are above the x-y-plane (we assume here a
 * positive solution, if it has negative values, you might want to decrease the
 * lower bound). Furthermore, it only takes the data points with z-values (#$3#)
 * equal to 0.5, i.e. a cut through the domain at #z=0.5#. For the data points
 * on this plane, the data values of the first data set (#$4#) are raised in
 * z-direction above the x-y-plane; all other points are denoted the value
 * #-1# instead of the value of the data vector and are not plotted due to
 * the lower bound in z plotting direction, given in the third pair of brackets.
 *
 * Of course, more complex cuts than #$3==.5# are possible, including nonlinear
 * ones. Note however, that only those points which are actually on the
 * cut-surface are plotted.
 *
 *
 * \subsection{POVRAY format}
 *
 * No information presently available.
 *
 *
 * \subsection{EPS format}
 *
 * To be filled in.
 * precision=5; viewpoint=gnuplot default; no border
 * shade or not; grid or not; grid shaded; data vector; memory; color=one per cell
 *
 * \subsection{GMV format}
 *
 * The #write_gmv# function and the #write# function through the #gmv# parameter
 * write the data in a format understood by the GMV (general mesh viewer)
 * program. This program is able to generate 2d and 3d plots of almost
 * arbitrarily many data sets, along with shading, cuts through data sets and
 * many other nifty features.
 *
 * Data is written in the following format: nodes are considered the points
 * of the patches. In spatial dimensions less than three, zeroes are
 * inserted for the missing coordinates. The data vectors are written as
 * node or cell data, where for the first the data space is interpolated to
 * (bi-,tri-)linear elements.
 *
 *
 * \section{Output parameters}
 *
 * All functions take a parameter which is a structure of type #XFlags#, where
 * #X# is the name of the output format. To find out what flags are presently
 * supported, read the documentation of the different structures.
 *
 * Note that usually the output formats used for scientific visualization
 * programs have no or very few parameters (apart from some compatibility flags)
 * because there the actual appearance of output is determined using the
 * visualization program and the files produced by this class store more or less
 * only raw data.
 *
 * The direct output formats, like Postscript or Povray need to be given a lot
 * more parameters, though, since there the output file has to contain all
 * details of the viewpoint, light source, etc.
 *
 *
 * @author Wolfgang Bangerth 1999; EPS output based on an earlier implementation by Stefan Nauber for the old DataOut class
 */
class DataOutBase 
{
  public:
				     /**
				      * Data structure describing a
				      * patch of data in #dim# space
				      * dimensions. See the general
				      * documentation of the
				      * #DataOutBase# class for more
				      * information on its contents
				      * and purposes.
				      *
				      * @author Wolfgang Bangerth */
    template <int dim>
    struct Patch
    {
					 /**
					  * Corner points of a patch.
					  * Inner points are computed by
					  * a multilinear transform of
					  * the unit cell to the cell
					  * specified by these corner
					  * points. The order of points
					  * is the same as for cells
					  * in the triangulation.
					  */
	Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
	
					 /**
					  * Number of subdivisions with
					  * which this patch is to be
					  * written. #1# means no
					  * subdivision, #2# means
					  * bisection, #3# trisection,
					  * etc.
					  */
	unsigned int n_subdivisions;
	
					 /**
					  * Data vectors. The format is
					  * as follows:
					  * #data(i,.)# denotes the data
					  * belonging to the #i#th data
					  * vector. #data.n()#
					  * therefore equals the number
					  * of output points; this
					  * number is #(subdivisions+1)^{dim}#.
					  * #data.m()# equals the number of
					  * data vectors.
					  *
					  * Within each column,
					  * #data(.,j)# are the data
					  * values at the output point #j#,
					  * where #j# runs first over the
					  * last direction, then over the second
					  * last one etc, just as if it was
					  * organized as an array
					  * #double[x][y][z]#. This is also
					  * the order of points as provided
					  * by the #QIterated# class when used
					  * with the #QTrapez# class as subquadrature.
					  * Note that if #subdivisions==1#,
					  * the elements of #data[i]# correspond
					  * to vertices #(0,1)# in 1d,
					  * #(0, 3, 1, 2)# in 2d, and
					  * #(0, 4, 3, 7, 1, 5, 2, 6)# in 3d.
					  *
					  * Since the number of data vectors
					  * is usually the same for all
					  * patches to be printed, #data.size()#
					  * should yield the same value for all
					  * patches provided.
					  */
	FullMatrix<double> data;
	
					 /**
					  * Default constructor.
					  */
	Patch () :
			n_subdivisions(0) {};
    };


				     /**
				      * Flags describing the details of
				      * output in UCD format.
				      */
    struct UcdFlags 
    {
					 /**
					  * Write a comment at the beginning
					  * of the file stating the date of
					  * creation and some other data.
					  * While this is supported by the
					  * UCD format (and the AVS program),
					  * some other programs get confused
					  * by this, so you can switch it off
					  * this way.
					  *
					  * Default: #true#.
					  */
	bool write_preamble;
	
					 /**
					  * Constructor.
					  */
	UcdFlags (const bool write_preamble = true);
    };

				     /**
				      * Flags describing the details of
				      * output in Gnuplot format. At
				      * present no flags are implemented.
				      */
    struct GnuplotFlags 
    {
      private:
					 /**
					  * Dummy entry to suppress compiler
					  * warnings when copying an empty
					  * structure. Remove this member
					  * when adding the first flag to
					  * this structure (and remove the
					  * #private# as well).
					  */
	int dummy;
    };

    				     /**
				      * Flags describing the details of
				      * output in Povray format. At
				      * present no flags are implemented.
				      */
    struct PovrayFlags 
    {
      private:
					 /**
					  * Dummy entry to suppress compiler
					  * warnings when copying an empty
					  * structure. Remove this member
					  * when adding the first flag to
					  * this structure (and remove the
					  * #private# as well).
					  */
	int dummy;
    };

    
    				     /**
				      * Flags describing the details of
				      * output in encapsulated postscript
				      * format.
				      */
    struct EpsFlags 
    {
					 /**
					  * This denotes the number of the
					  * data vector which shall be used
					  * for generating the height
					  * information. By default, the
					  * first data vector is taken,
					  * i.e. #height_vector==0#, if
					  * there is any data vector. If there
					  * is no data vector, no height
					  * information is generated.
					  */
	unsigned int height_vector;

					 /**
					  * Number of the vector which is
					  * to be taken to colorize cells.
					  * The same applies as for
					  * #height_vector#.
					  */
	unsigned int color_vector;
	
					 /**
					  * Enum denoting the possibilities
					  * whether the scaling should be done
					  * such that the given #size# equals
					  * the width or the height of
					  * the resulting picture.
					  */
	enum SizeType {
	      width, height
	};

					 /**
					  * See above. Default is #width#.
					  */
	SizeType size_type;
	
					 /**
					  * Width or height of the output
					  * as given in postscript units
					  * This usually is given by the
					  * strange unit 1/72 inch. Whether
					  * this is height or width is
					  * specified by the flag
					  * #size_type#.
					  *
					  * Default is 300.
					  */
	unsigned int size;

					 /**
					  * Width of a line in postscript
					  * units. Default is 0.5.
					  */
	double line_width;
	
					 /**
					  * Angle of the line origin-viewer
					  * against the z-axis in degrees.
					  *
					  * Default is the Gnuplot-default
					  * of 60.
					  */
	double azimut_angle;

					 /**
					  * Angle by which the viewers
					  * position projected onto the
					  * x-y-plane is rotated around
					  * the z-axis, in positive sense
					  * when viewed from above. The
					  * unit are degrees, and zero
					  * equals a position above or below
					  * the negative y-axis.
					  *
					  * Default is the Gnuplot-default
					  * of 30.
					  */
	double turn_angle;

					 /**
					  * Factor by which the z-axis is to
					  * be stretched as compared to the
					  * x- and y-axes. This is to compensate
					  * for the different sizes that
					  * coordinate and solution values may
					  * have and to prevent that the plot
					  * looks to much out-of-place (no
					  * elevation at all if solution values
					  * are much smaller than coordinate
					  * values, or the common "extremely
					  * mountainous area" in the opposite
					  * case.
					  *
					  * Default is #1.0#.
					  */
	double z_scaling;

					 /**
					  * Flag the determines whether the
					  * lines bounding the cells (or the
					  * parts of each patch) are to be
					  * plotted.
					  *
					  * Default: #true#.
					  */
	bool   draw_mesh;

					 /**
					  * Flag whether to fill the regions
					  * between the lines bounding the cells
					  * or not. If not, no hidden line removal
					  * is performed, which in this crude
					  * implementation is done through
					  * writing the cells in a back-to-front
					  * order, thereby hiding the cells in
					  * the background by cells in the
					  * foreground.
					  *
					  * If this flag is #false# and #draw_mesh#
					  * is #false# as well, nothing will be
					  * printed.
					  *
					  * If this flag is #true#, then the cells
					  * will be drawn either colored by one
					  * of the data sets (if #shade_cells# is
					  * #true#), or pure white (if
					  * #shade_cells# is false or if there are
					  * no data sets).
					  *
					  * Default is #true#.
					  */
	bool   draw_cells;

					 /**
					  * Flag to determine whether the cells
					  * shall be colorized by one the data
					  * set denoted by #color_vector#, or
					  * simply be painted in white. This
					  * flag only makes sense if
					  * #draw_cells==true#. Colorization is
					  * done through the #color_function#.
					  *
					  * Default is #true#.
					  */
	bool   shade_cells;

					 /**
					  * Structure keeping the three color
					  * values in the RGB system.
					  */
	struct RgbValues { float red, green, blue; };

					 /**
					  * Definition of a function pointer
					  * type taking a value and returning
					  * a triple of color values in RGB
					  * values.
					  *
					  * Besides the actual value by which
					  * the color is to be computed, min
					  * and max values of the data to
					  * be colorized are given as well.
					  */
	typedef RgbValues (*ColorFunction) (const double value,
					    const double min_value,
					    const double max_value);

					 /**
					  * This is a pointer to the function
					  * which is used to colorize the cells.
					  * By default, it points to the
					  * static function #default_color_function#
					  * which is a member of this class.
					  */
	ColorFunction color_function;


					 /**
					  * Default colorization function. This
					  * one does what one usually wants:
					  * it shifts colors from black (lowest
					  * value) through blue, green and red
					  * to white (highest value). For the
					  * exact defition of the color scale
					  * refer to the impementation.
					  *
					  * This function was originally written
					  * by Stefan Nauber.
					  */
	static RgbValues default_color_function (const double value,
						 const double min_value,
						 const double max_value);
	
					 /**
					  * Constructor.
					  */
	EpsFlags (const unsigned int  height_vector = 0,
		  const unsigned int  color_vector  = 0,
		  const SizeType      size_type     = width,
		  const unsigned int  size          = 300,
		  const double        line_width    = 0.5,
		  const double        azimut_angle  = 60,
		  const double        turn_angle    = 30,
		  const double        z_scaling     = 1.0,
		  const bool          draw_mesh     = true,
		  const bool          draw_cells    = true,
		  const bool          shade_cells   = true,
		  const ColorFunction color_function= &default_color_function);
    };

    				     /**
				      * Flags describing the details of
				      * output in gmv format. At
				      * present no flags are implemented.
				      */
    struct GmvFlags 
    {
      private:
					 /**
					  * Dummy entry to suppress compiler
					  * warnings when copying an empty
					  * structure. Remove this member
					  * when adding the first flag to
					  * this structure (and remove the
					  * #private# as well).
					  */
	int dummy;
    };

				     /**
				      * Write the given list of patches
				      * to the output stream in ucd
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim>
    static void write_ucd (const vector<Patch<dim> > &patches,
			   const vector<string>      &data_names,
			   const UcdFlags            &flags,
			   ostream                   &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in gnuplot
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim>
    static void write_gnuplot (const vector<Patch<dim> > &patches,
			       const vector<string>      &data_names,
			       const GnuplotFlags        &flags,
			       ostream                   &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in povray
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim>
    static void write_povray (const vector<Patch<dim> > &patches,
			      const vector<string>      &data_names,
			      const PovrayFlags         &flags,
			      ostream                   &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in eps
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim>
    static void write_eps (const vector<Patch<dim> > &patches,
			   const vector<string>      &data_names,
			   const EpsFlags            &flags,
			   ostream                   &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in gmv
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim>
    static void write_gmv (const vector<Patch<dim> > &patches,
			   const vector<string>      &data_names,
			   const GmvFlags            &flags,
			   ostream                   &out);

    
				     /**
				      * Exception
				      */
    DeclException2 (ExcUnexpectedNumberOfDatasets,
		    int, int,
		    << "The number of data sets on this patch is " << arg1
		    << ", but we expected " << arg2);

				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidDatasetSize,
		    int, int,
		    << "The number of points in this data set is " << arg1
		    << ", but we expected " << arg2 << " in each space direction.");
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidVectorNumber,
		    int, int,
		    << "The number " << arg1 << " of the vector to be used for "
		    << "this information is invalid, since there are only "
		    << arg2 << " data sets.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);

  private:
				     /**
				      * Class holding the data of one
				      * cell of a patch in two space
				      * dimensions for output. It is
				      * the projection of a cell in
				      * three-dimensional space (two
				      * coordinates, one height value)
				      * to the direction of sight.
				      */
    class EpsCell2d {
      public:
	
					 /**
					  * Vector of vertices of this cell.
					  */
	Point<2> vertices[4];
	
					 /**
					  * Data value from which the actual
					  * colors will be computed by
					  * the colorization function stated
					  * in the #EpsFlags# class.
					  */
	float color_value;
	
					 /**
					  * Depth into the picture, which
					  * is defined as the distance from
					  * an observer at an the origin in
					  * direction of the line of sight.
					  */
	float depth;
	
					 /**
					  * Comparison operator for
					  * sorting.
					  */
	bool operator < (const EpsCell2d &) const;
    };

};

	




/**
 * This class is the interface to the #DataOutBase# class, as already its name
 * might suggest. It does not offer much functionality apart from a way to
 * access the implemented formats and a way to dynamically dispatch what output
 * format to chose.
 *
 * This class is thought as a base class to classes actually
 * generating data for output. It has two abstract virtual functions,
 * #get_patches# and #get_dataset_names# which are to produce the data
 * which is actually needed. These are the only functions that need to
 * be overloaded by a derived class.  In additional to that, it has a
 * function for each output format supported by the underlying base
 * class which gets the output data using these two virtual functions
 * and passes them to the raw output functions.
 *
 * The purpose of this class is mainly two-fold: to support storing flags
 * by which the output in the different output formats are controlled,
 * and means to work with output in a way where output format, flags and
 * other things are determined at run time. In addition to that it offers
 * the abstract interface to derived classes briefly discussed above.
 *
 *
 * \subsection{Output flags}
 *
 * The way we treat flags in this class is very similar to that used in
 * the #GridOut# class. For detailed information on the why's and how's,
 * as well as an example of programming, we refer to the documentation
 * of that class.
 *
 * In basics, this class stores a set of flags for each output format
 * supported by the underlying #DataOutBase# class. These are used
 * whenever one of the #write_*# functions is used. By default, the
 * values of these flags are set to reasonable start-ups, but in case
 * you want to change them, you can create a structure holding the flags
 * for one of the output formats and set it using the #set_flags# functions
 * of this class to determine all future output the object might produce
 * by that output format.
 *
 * For information on what parameters are supported by different output
 * functions, please see the documentation of the #DataOutBase# class and
 * its member classes.
 *
 *
 * \subsection{Run time selection of formats}
 *
 * This class, much like the #GridOut# class, has a set of functions
 * providing a list of supported output formats, an #enum# denoting
 * all these and a function to parse a string and return the respective
 * #enum# value if it is a valid output format's name. Finally, there
 * is a function #write#, which takes a value of this #enum# and
 * dispatches to one of the actual #write_*# functions depending on
 * the output format selected by this value. 
 *
 * The functions offering the different output format names are,
 * respectively, #default_suffix#, #parse_output_format#, and
 * #get_output_format_names#. They make the selection of ouput formats
 * in parameter files much easier, and especially independent of
 * the formats presently implemented. User programs need therefore not
 * be changed whenever a new format is implemented.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class DataOutInterface : private DataOutBase
{
  public:
				     /**
				      * Provide a data type specifying the
				      * presently supported output formats.
				      */
    enum OutputFormat { ucd, gnuplot, povray, eps, gmv };

				     /**
				      * Obtain data through the #get_patches#
				      * function and write it to #out# in
				      * UCD format.
				      */
    void write_ucd (ostream &out) const;

				     /**
				      * Obtain data through the #get_patches#
				      * function and write it to #out# in
				      * GNUPLOT format.
				      */
    void write_gnuplot (ostream &out) const;

    				     /**
				      * Obtain data through the #get_patches#
				      * function and write it to #out# in
				      * POVRAY format.
				      */
    void write_povray (ostream &out) const;

				     /**
				      * Obtain data through the #get_patches#
				      * function and write it to #out# in
				      * EPS format.
				      */
    void write_eps (ostream &out) const;

    				     /**
				      * Obtain data through the #get_patches#
				      * function and write it to #out# in
				      * GMV format.
				      */
    void write_gmv (ostream &out) const;

				     /**
				      * Write data and grid to #out# according
				      * to the given data format. This function
				      * simply calls the appropriate
				      * #write_*# function.
				      */
    void write (ostream &out, const OutputFormat output_format) const;
    
				     /**
				      * Set the flags to be used for output
				      * in UCD format.
				      */
    void set_flags (const UcdFlags &ucd_flags);

    				     /**
				      * Set the flags to be used for output
				      * in GNUPLOT format.
				      */
    void set_flags (const GnuplotFlags &gnuplot_flags);

    				     /**
				      * Set the flags to be used for output
				      * in POVRAY format.
				      */
    void set_flags (const PovrayFlags &povray_flags);

    				     /**
				      * Set the flags to be used for output
				      * in 1d EPS output.
				      */
    void set_flags (const EpsFlags &eps_flags);

    				     /**
				      * Set the flags to be used for output
				      * in GMV format.
				      */
    void set_flags (const GmvFlags &gmv_flags);

    
				     /**
				      * Provide a function which tells us which
				      * suffix with a given output format
				      * usually has. At present the following
				      * formats are defined:
				      * \begin{itemize}
				      * \item #ucd#: #.inp#
				      * \item #gnuplot#: #.gnuplot#
				      * \item #povray#: #.pov#
				      * \item #eps#: #.eps#
				      * \item #gmv#: #.gmv#.
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
    
  protected:
				     /**
				      * This is the abstract function through
				      * which derived classes propagate
				      * preprocessed data in the form of
				      * #Patch# structures (declared in
				      * the base class #DataOutBase#) to
				      * the actual output function. You
				      * need to overload this function to
				      * allow the output functions to
				      * know what they shall print.
				      */
    virtual const vector<DataOutBase::Patch<dim> > & get_patches () const = 0;

				     /**
				      * Abstract virtual function through
				      * which the names of data sets are
				      * obtained by the output functions
				      * of the base class.
				      */
    virtual vector<string> get_dataset_names () const = 0;

  private:
				     /**
				      * Flags to be used upon output of UCD
				      * data. Can be changed by using the
				      * #set_flags# function.
				      */
    UcdFlags     ucd_flags;

				     /**
				      * Flags to be used upon output of GNUPLOT
				      * data. Can be changed by using the
				      * #set_flags# function.
				      */
    GnuplotFlags gnuplot_flags;

    				     /**
				      * Flags to be used upon output of POVRAY
				      * data. Can be changed by using the
				      * #set_flags# function.
				      */
    PovrayFlags povray_flags;

				     /**
				      * Flags to be used upon output of EPS
				      * data in one space dimension. Can be
				      * changed by using the #set_flags#
				      * function.
				      */
    EpsFlags     eps_flags;    

				     /**
				      * Flags to be used upon output of gmv
				      * data in one space dimension. Can be
				      * changed by using the #set_flags#
				      * function.
				      */
    GmvFlags     gmv_flags;    
};




/*----------------------------   data_out_base.h     ---------------------------*/
/* end of #ifndef __data_out_base_H */
#endif
/*----------------------------   data_out_base.h     ---------------------------*/
