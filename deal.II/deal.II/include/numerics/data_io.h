//----------------------------  data_io.h  ---------------------------
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
//----------------------------  data_io.h  ---------------------------
#ifndef __deal2__data_io_h
#define __deal2__data_io_h


#include <base/exceptions.h>
#include <vector>
#include <string>

template <typename number> class Vector;
template <int dim> class DoFHandler;
class EpsOutputData;

/**
 * This class is deprecated. Use the #DataOut# class instead.
 *
 * This class implements an output mechanism for grid and simulation data
 * in several formats.
 * At present it supports output in UCD (unstructured cell data) and
 * GNUPLOT format. Partly supported are POVRAY and encapsulated postscript.
 * POVRAY allows for only one data set and does not support cell data.
 * Encapsulated Postscript supports Q1-Elements only.
 *
 * It allows the user to attach a degree of freedom handler object
 * (#DoFHandler#) which also gives access to the geometry data of the
 * underlying triangulation and to add data vectors of which the values
 * are to be written.
 *
 *
 * \subsection{Limitations}
 * 
 * At present, no grouping of components to vectors is implemented, i.e.
 * you can only write each component independent of the others. Also, it
 * is not possible to output calculations which were performed on elements
 * with more or less than one degree of freedom per vertex.
 *
 * 
 * \subsection{UCD format}
 *
 * The UCD format is described in the AVS developer's guide. Due to
 * limitations in the present format, only node based data can be output,
 * so higher order elements are only written with their node values, no
 * interior or line values are used. No use is made of the possibility
 * to give model data since these are not supported by all
 * UCD aware programs. You may give cell data, but that is not supported by
 * all programs.
 *
 * The ASCII UCD format is used. In future versions, a binary version may
 * follow up.
 *
 * Note that to enumerate the vertices, not the vertex index is used but
 * the index of the degree of freedom located on this vertex. This makes
 * the mapping between the vertices and the entries in the data vectors
 * much easier.
 *
 *
 * \subsection{GNUPLOT draft format}
 *
 * The GNUPLOT format is not able to handle data on unstructured grids
 * directly. Directly would mean that you only give the vertices and
 * the solution values thereon and the program constructs its own grid
 * to represent the data. This is only possible for a structured tensor
 * product grid in two dimensions.
 *
 * In one dimension, the format is obviously #x v1 v2 ...#, where #x#
 * is the coordinate value of a grid point, while the #vi# are the
 * vector elements referring to the present node. Within GNUPLOT,
 * call #plot "filename" using 1:x#. #x# denotes the number of the data set you
 * want to see plus one. For example #using 1:4# would mean to plot the
 * third data vector.
 *
 * For more than one dimension, the #DataOut_Old<dim>::write_gnuplot()# somehow
 * duplicates the functionality of the #Triangulation<dim>::print_gnuplot()#
 * functions. These, however, offer more functionality in some respect.
 * The grid is represented as a sequence of lines, where each cell is
 * a sequence of five vertices (the first one is appended at the end to close
 * the contour of the cell) with the data appended after each vertex. Each cell
 * is therefore a sequence of five lines #x y v1 v2 ...# forming together
 * the bounding line of this cell. After each cell, two newlines are inserted
 * to prevent GNUPLOT from joining the lines bounding two cells.
 *
 * To view the results in two dimensions, use #set data style lines#
 * within gnuplot and call #plot "filename"# to see the grid. Use
 * #set parametric# and #splot "filename" using 1:2:x# to get a 3d surface
 * plot of the (#x-2#)th data set. For example, using #x=4# would mean to
 * plot the second data set.
 *
 * This format is somewhat restricted to continuous data and to finite elements
 * of first order only. The reason for the first restriction is that it takes
 * nodal values and can therefore only work if a finite element has degrees of
 * freedom in the vertices of each cell. This is not the case for discontinuous
 * elements. The second restriction is only a problem of quality of output: you
 * actually can print quadratic or higher order elements using this style, but
 * you will only seem the contour of each cell with the bounding lines of each
 * cell being straight lines. You won't see the structure of the solution in the
 * interior of a cell nor on the lines bounding it.
 *
 *
 * \subsection{GNUPLOT 'quality' format}
 *
 * To remedy the abovementioned problems, the GNUPLOT 'quality' format was
 * introduced. As noted above, GNUPLOT can handle tensor grids and displays
 * them quite nicely, including hidden line removal. It can also handle more
 * than one tensor grid, then called a patch. The idea now is to describe
 * each cell as a patch of $N\times N$ points (in 2D). For linear elements,
 * $2\times 2$ is sufficient, but for higher elements you will want to use
 * a higher number. The #write_gnuplot# function writes the data in this format
 * and the argument it takes is the number of subintervals it divides each cell
 * into, i.e. $N-1$.
 *
 * This output routine also addresses the problem introduced with discontinuous
 * elements, since it takes its data locally from each cell and displays it
 * as a patch separately; it therefore does not use continuity and GNUPLOT will
 * also plot discontinuities, if they are there.
 *
 *
 * \subsection{POVRAY mesh format}
 *
 * POVRAY is a ray tracing tool in the public domain and poduces high quality
 * images of three dimensional scenes. Since this tool can handle only one
 * scene at a time, only the first data vector is output, but this may change
 * in future versions of this library. For a reference of the data format
 * and the usage of POVRAY see the extensive manual of that program.
 *
 * The POVRAY output format is presently only supported for two spatial
 * dimensions. It displays every quadrilateral as two triangles, which clearly
 * is not an optimal solution, but necessary since POVRAY presently does not
 * support bilinear quadrilaterals and using polygons as vastly suboptimal in
 * term of memory and speed to the triangle mesh supported by POVRAY.
 *
 *
 * \subsection{Encapsulated Postscript format}
 *
 * There is a function for generating encapsulated Postscript
 * (EPS) without the need for another graphics tool. This functionality
 * is provided through #write_eps# and #EpsOutputData#. 
 * #write_eps# can use one data vector for height information and one
 * cell vector for shading information. Control is done through an
 * object of class #EpsOutputData# as follows or by default.
 *
 * Vectors are added as usual by #add_data_vector#. Then one has to
 * decide, wether to produce a 2D or 3D plot. This is done by setting
 * #height_info# to 
 * \begin{description} 
 *   \item[NoHeight] for 2D-Output (or Top-View thats the same by no
 *     turning is done) or to
 *   \item[HeightVector] for 3D-Output. You have to attach a
 *     #dof_data_vector# to actually get 3D. If you don't then output
 *     will be generated in 2D.
 *   \item[DefaultHeight] is 3D if there is a #dof_data_vector# and 2D if
 *     none is present.
 * \end{description}
 * For 3D-Output one has to set #azimuth# and #elevation# for the
 * angle of view and #height_vector# to the number of the #dof_data#
 * vector that provides the height information to be used. The default
 * values are analogous to GNUPlot for azimuth and elevation and
 * vector 0 for the height vector.
 *
 * The cells can be shaded in four different modes, controlled by the
 * attribute #cell_shading#:
 * \begin{enumerate}
 *   \item[NoShading] provides transparent shading.
 *   \item[ShadingVector] uses a cell vector to do shading. The number
 *     of the cell vector to be uses is provided in #cell_vector#. To
 *     scale the cell vector there is the method #color#. It is called
 *     with the actual value of the cell, the maximum and the minimum
 *     value of a cell in the cell vector. It returns three values for
 *     red, green and blue. If there no #cell_data# vector than there is
 *     transparent shading.
 *   \item[LightShaded] just shades the plot. This is controlled by
 *     the vector #light# which stores the direction of the light 
 *     beams. This is done only if there is height information.
 *   \item[DefaultShading] is controlled by presence of different
 *     vectors. If there no height information then do no
 *     shading. Otherwise if there is #cell_data# use this for shading.
 *     Otherwise do light shading.
 * \end{enumerate}
 *
 * Finnaly one can choose to mark the cell boundaries by setting
 * #cell_boundary_shading#. It can take one of four values:
 * \begin{itemize}
 *   \item NoBoundary for no cell boundaries,
 *   \item DefaultBoundary or
 *   \item BlackBoundary for black cell boundaries,
 *   \item WhiteBoundary for white cell boundaries, 
 * \end{itemize}
 *
 * Another interesting feature is that you can write multiple
 * eps-pictures to one file by just doing several invocations of
 * #write_eps#. One than can switch between the different graphics
 * with the #>># Button in GhostView for example.
 *
 *
 * \subsection{GMV format}
 *
 * The #write_gmv# function and the #write# function through the #gmv# parameter
 * write the data in a format understood by the GMV (general mesh viewer)
 * program. This program is able to generate 2d and 3d plots of almost
 * arbitrarily many data sets, along with shading, cuts through data sets and
 * many other nifty features.
 *
 * Data is written in the following format: nodes are considered the vertices
 * of the triangulation. In spatial dimensions less than three, zeroes are
 * inserted for the missing coordinates. The data vectors are written as
 * node or cell data, where for the first the data space is interpolated to
 * (bi-,tri-)linear elements.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, Stefan Nauber, 1998, 1999
 */
template <int dim>  
class DataOut_Old {
  public:
				     /**
				      * Provide a data type specifying the
				      * presently supported output formats.
				      */
    enum OutputFormat { ucd, gnuplot, gnuplot_draft, povray_mesh, eps, gmv };
    
				     /**
				      * Constructor
				      */
    DataOut_Old ();
    
				     /**
				      * Designate a dof handler to be used
				      * to extract geometry data and the
				      * mapping between nodes and node values.
				      */
    void attach_dof_handler (const DoFHandler<dim> &);

				     /**
				      * Add a data vector together with its
				      * name and the physical unit
				      * (for example meter, kelvin, etc). By
				      * default, "<dimensionless>" is assumed
				      * for the units.
				      *
				      * A pointer to the vector is stored, so
				      * you have to make sure the vector
				      * exists at that address at least as
				      * long as you call the
				      * #write_*# functions.
				      *
				      * It is assumed that the vector has the
				      * same number of components as there are
				      * degrees of freedom in the dof handler,
				      * in which case it is assumed to be a
				      * vector storing nodal data; or the size
				      * may be the number of active cells on
				      * the present grid, in which case it is
				      * assumed to be a cell data vector.
				      *
				      * The name and unit of a data vector shall
				      * only contain characters which are
				      * letters, underscore and a few other
				      * ones. Refer to the #ExcInvalidCharacter#
				      * exception declared in this class to
				      * see which characters are valid and which
				      * are not.
				      */
    void add_data_vector (const Vector<double> &data,
			  const string         &name,
			  const string         &units="<dimensionless>");

				     /**
				      * Release the pointers to the data
				      * vectors. You have to set all data
				      * entries again using the
				      * #add_data_vector# function. The pointer
				      * to the dof handler remains stored,
				      * however.
				      */
    void clear_data_vectors ();

				     /**
				      * Write the stored data to the given
				      * stream in UCD data. You may have
				      * written any comment to that stream
				      * before calling this function. Comments
				      * start with the \# character in the
				      * first column of a line and may only
				      * appear at the beginning of a file,
				      * without non-comment lines inbetween.
				      */
    void write_ucd (ostream &out) const;

				     /**
				      * Write data and grid in GNUPLOT format.
				      * Write a patch for each grid cell using
				      * #accuracy# subintervals per dimension.
				      */
    void write_gnuplot (ostream &out, const unsigned int accuracy=1) const;

				     /**
				      * Write data and grid in GNUPLOT format.
				      * Only the edges between cells are written.
				      */
    void write_gnuplot_draft (ostream &out) const;

				     /**
				      * Write data of first vector and grid
				      * in POVRAY format. Further data vectors
				      * as well as cell data is ignored.
				      */
    void write_povray_mesh (ostream &out) const;

				   /**
				    * Write data of first vector and grid in
				    * Encapsulated postscript format.
				    * Further data vectors
				    * as well as cell data is ignored.
				    *
				    * If no special output data is given, a
				    * default constructed object is used.
				    */
    void write_eps (ostream &out,
		    const EpsOutputData &eod = EpsOutputData()) const;

    				   /**
				    * Write data and grid in GMV format.
				    */
    void write_gmv (ostream &out) const;

				     /**
				      * Write data and grid to #out# according
				      * to the given data format. This function
				      * simply calls the appropriate
				      * #write_*# function.
				      *
				      * In case of gnuplot output, the
				      * standard accuracy of 1 is
				      * chosen.
				      */
    void write (ostream &out, const OutputFormat output_format) const;
    
				     /**
				      * Provide a function which tells us which
				      * suffix with a given output format
				      * usually has. At present the following
				      * formats are defined:
				      * \begin{itemize}
				      * \item #ucd#: #.inp#
				      * \item #gnuplot# and #gnuplot_draft#:
				      *    #.gnuplot#
				      * \item #povray_mesh#: #.pov#
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
    DeclException0 (ExcIncorrectDofsPerVertex);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoDoFHandlerSelected);
				     /**
				      * Exception
				      */
    DeclException3 (ExcInvalidVectorSize,
		    int, int, int,
		    << "The vector has size " << arg1
		    << " but the DoFHandler objects says there are " << arg2
		    << " degrees of freedom and there are " << arg3
		    << " active cells.");
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidCharacter,
		    string,
		    << "Please use only the characters [a-zA-Z0-9_<>()] for" << endl
		    << "description strings since AVS will only accept these." << endl
		    << "The string you gave was <" << arg1 << ">.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidState);
    
  private:

				     /**
				      * Declare an entry in the list of
				      * data elements.
				      */
    struct DataEntry {
					 /**
					  * Default constructor, constructs
					  * invalid object, but is needed for
					  * the STL classes.
					  */
	DataEntry ();
			
					 /**
					  * Constructor
					  */
	DataEntry (const Vector<double> *data, const string name, const string units);
	
					 /**
					  * Pointer to the data vector.
					  */
	const Vector<double> *data;
	
					 /**
					  * Name of this component.
					  */
	string   name;
	
					 /**
					  * Physical unit name of this
					  * component.
					  */
	string   units;
    };

				     /**
				      * Pointer to the dof handler object.
				      */
    const DoFHandler<dim>   *dofs;

				     /**
				      * List of data elements with vectors of
				      * values for each degree of freedom.
				      */
    vector<DataEntry>  dof_data;

				     /**
				      * List of data elements with vectors of
				      * values for each cell.
				      */
    vector<DataEntry>  cell_data;


/**
				      * Return the number of faces in the
				      * triangulation which have a boundary
				      * indicator not equal to zero. Only
				      * these faces are explicitely printed
				      * in the #write_*# functions;
				      * all faces with indicator 255 are
				      * interior ones and an indicator with
				      * value zero for faces at the boundary
				      * are considered default.
				      *
				      * This function always returns an empty
				      * list in one dimension.
				      *
				      * The reason for this function is the
				      * same as for #write_ucd_faces#. See
				      * there for more information.
				      */
    unsigned int n_boundary_faces () const;

				     /**
				      * Write the grid information about
				      * faces to #out#. Only those faces
				      * are printed which are on the boundary
				      * and which have a boundary indicator
				      * not equal to zero, since the latter
				      * is the default for boundary faces.
				      *
				      * Since cells and faces are continuously
				      * numbered, the #starting_index# for
				      * the numbering of the faces is passed
				      * also.
				      *
				      * This function unfortunately can not
				      * be included in the regular #write_ucd#
				      * function, since it needs special
				      * treatment for the case #dim==1#, in
				      * which case the face iterators are
				      * #void*#'s and lack the member functions
				      * which are called. We would not actually
				      * call these functions, but the compiler
				      * would complain anyway when compiling
				      * the function for #dim==1#. Bad luck.
				      */
    void write_ucd_faces (ostream &out,
			  const unsigned int starting_index) const;

				     /**
				      * Class of data for vertices
				      * for eps output.
				      */
    class EpsVertexData{
      public:
				       /**
					* Coordinates
					*/
        double x;
	double y;
	double z;
      
				       /**
					* Default Constructor;
					* needed for STL.
					*/
        EpsVertexData()
	  {};
      
				       /**
					* Constructor that takes three
					* coordinate values as argument.
					*/
        EpsVertexData(double a, double b, double c):x(a),y(b),z(c)
	  {};
      
				       /**
					* Transformation method
					*/
	void turn(double azi, double ele);
    };
    
                                     /** 
				      * Class of cell data for output.
				      * For eps output some calculations have
				      * to be done between transformation and
				      * projection and output. Because of this all
				      * output data is put into a STL multiset.
				      */
    class EpsCellData{
      public:
	
					 /**
					  * STL-vector of vertices
					  */
	EpsVertexData vertices[4];
	
					 /**
					  * Color values
					  */
	float red;
	float green;
	float blue;

					 /**
					  * Transformation method
					  */
	void turn(double azi, double ele);
	
					 /**
					  * Comparison operator for
					  * sorting
					  */
	bool operator < (const EpsCellData &) const;
    };
};


/**
 * Structure for the control of encapsulated postscript output. See
 * general documentation of class #DataOut_Old# for description.
 *
 * @author Stefan Nauber
 */
class EpsOutputData{
  public:
				     /**
				      * Types of height info
				      */
    enum HeightInfo {
	  DefaultHeight, NoHeight, HeightVector
    };
    
				     /**
				      * Types of cell shading
				      */
    enum CellShading {
	  DefaultShading, NoShading, ShadingVector, LightShaded
    };

				     /**
				      * Types of cell boundary shading
				      */
    enum CellBoundaryShading {
	  DefaultBoundary, NoBoundary, BlackBoundary, WhiteBoundary
    };

				     /**
				      * Where height information comes from
				      */
    HeightInfo height_info;

				     /** 
				      * If and how cells are shaded
				      */
    CellShading cell_shading;

				     /**
				      * Color selection when shading cell
				      * boundaries.
				      */
    CellBoundaryShading cell_boundary_shading;

				     /**
				      * Number of the vector which is to be
				      * used for the height information, within
				      * the list of DoF data vectors.
				      */
    unsigned int height_vector;

				     /**
				      * Number of the vector which is to be
				      * used for the cell shading values, within
				      * the list of cell data vectors.
				      */
    unsigned int cell_vector;
    
                                     /**
				      * Azimuth of the spectators position.
				      * This defines the position on the
				      * x-y-plane and is an angle of rotation
				      * around the z-axis. We define that
				      * if the azimuth angle is zero, this
				      * means that the observer is sitting on
				      * (or above or below) the positive y-axis
				      * and looks back to the origin (direction
				      * of view is always towards the origin).
				      * Positive angles denote rotation in
				      * clockwise sense (as viewed from the top),
				      * i.e. 90 degrees would be sitting on the
				      * positive x-axis, etc.
				      *
				      * Please note that the angle has to be
				      * given in radians rather than in degrees.
				      */
    float azimuth;

                                     /**
				      * Elevation of the spectators position.
				      * This is the angle that the line between
				      * the origin and the spectators position
				      * forms with the x-y-plane, measured
				      * upwards in direction towards the
				      * positive z-axis. 
				      *
				      * Please note that the angle has to be
				      * given in radians rather than in degrees.
				      */
    float elevation;

				     /**
				      * Direction of the light beams. 
				      */
    double light[3];
    
				     /**
				      * Default constructor. Sets height and
				      * shading flags to their default values,
				      * azimut and elevation to 0.2 each, and
				      * puts the light source at #(-1,-1,0)#.
				      */
    EpsOutputData();

				     /**
				      * Function returning a color value in
				      * rgb variables corresponding to the
				      * given value #x#. #x# is considered
				      * to be a values within the range
				      * #xmin...xmax#.
				      */
    void color (const float x,
		const float xmax,
		const float xmin, 
		float      &r,
		float      &g,
		float      &b) const;
};


#endif

