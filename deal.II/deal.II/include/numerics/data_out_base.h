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
 * \subsectin{Patches}
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
 * #(x0,y0) (x1,y0) (x2,y0) (x0,y1) (x1,y1) (x2,y1) (x0,y2) (x1,y2) (x2,y2)#,
 * i.e. the x-coordinate runs fastest, then the y-coordinate, then z (if there
 * are that many space directions).
 *
 * The #Patch# class takes a template parameter denoting the space dimension
 * in which this patch operates.
 *
 *
 * @author Wolfgang Bangerth 1999
 */
class DataOutBase 
{
  public:

				     /**
				      * Data structure describing a
				      * patch of data in #dim# space
				      * dimensions. See the general
				      * documentation for more information
				      * on its contents and purposes.
				      */
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
	const Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];

					 /**
					  * Number of subdivisions with
					  * which this patch is to be
					  * written. #1# means no
					  * subdivision, #2# means
					  * bisection, #3# trisection,
					  * etc.
					  */
	const unsigned int n_subdivisions;

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
					  * direction spanned by the first
					  * two points of the #corners#
					  * array, then over the direction
					  * spanned by points 0 and 3 and
					  * finally that spanned by points
					  * 0 and 4 (for 3d, for lower
					  * dimensions, this row is
					  * truncated, of course).
					  * Note that if #subdivisions==1#,
					  * the elements of #data[i]# correspond
					  * to vertices #(0,1)# in 1d,
					  * #(0, 1, 3, 2)# in 2d, and
					  * #(0, 1, 3, 2, 4, 5, 7, 6)# in 3d.
					  *
					  * Since the number of data vectors
					  * is usually the same for all
					  * patches to be printed, #data.size()#
					  * should yield the same value for all
					  * patches provided.
					  */
	const FullMatrix<double> data;

					 /**
					  * Default constructor.
					  */
	Patch ();
    };


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
};

	


/*----------------------------   data_out_base.h     ---------------------------*/
/* end of #ifndef __data_out_base_H */
#endif
/*----------------------------   data_out_base.h     ---------------------------*/
