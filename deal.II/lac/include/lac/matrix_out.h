//----------------------------  matrix_out.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  matrix_out.h  ---------------------------
#ifndef __deal2__matrix_out_h
#define __deal2__matrix_out_h

#include <base/data_out_base.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>


/**
 * Output a matrix in graphical form using the generic format
 * independent output routines of the base class. The matrix is
 * converted into a list of patches on a 2d domain where the height is
 * given by the elements of the matrix. The functions of the base
 * class can then write this "mountain representation" of the matrix
 * in a variety of graphical output formats. The coordinates of the
 * matrix output are that the columns run with increasing x-axis, as
 * usual, starting from zero, while the rows run into the negative
 * y-axis, also starting from zero. Note that due to some internal
 * restrictions, this class can only output one matrix at a time,
 * i.e. it can not take advantage of the multiple dataset capabilities
 * of the base class.
 *
 * A typical usage of this class would be as follows:
 * @begin{verbatim}
 *    FullMatrix<double> M;
 *    ...                // fill matrix M with some values
 *
 *                       // now write out M:
 *    MatrixOut matrix_out;
 *    std::ofstream out ("M.gnuplot");
 *    matrix_out.build_patches (M, "M");
 *    matrix_out.write_gnuplot (out);
 * @end{verbatim}
 * Of course, you can as well choose a different graphical output
 * format. Also, this class supports any matrix, not only of type
 * @ref{FullMatrix}, as long as it satisfies a number of requirements,
 * stated with the member functions of this class.
 *
 * @author Wolfgang Bangerth, 2001
 */
class MatrixOut : public DataOutInterface<2,2>
{
  public:
				     /**
				      * Destructor. Declared in order
				      * to make it virtual.
				      */
    virtual ~MatrixOut ();
    
				     /**
				      * Generate a list of patches
				      * from the given matrix and use
				      * the given string as the name
				      * of the data set upon writing
				      * to a file. Once patches have
				      * been built, you can use the
				      * functions of the base class to
				      * write the data into a files,
				      * using one of the supported
				      * output formats.
				      *
				      * Note that this function
				      * requires that we can extract
				      * elements of the matrix, which
				      * is done using the
				      * @p{get_element} function
				      * declared below. By adding
				      * specializations, you can
				      * extend this class to other
				      * matrix classes which are not
				      * presently
				      * supported. Furthermore, we
				      * need to be able to extract the
				      * size of the matrix, for which
				      * we assume that the matrix type
				      * offers member functions
				      * @p{m()} and @p{n()}, which
				      * return the number of rows and
				      * columns, respectively.
				      */
    template <class Matrix>
    void build_patches (const Matrix      &matrix,
			const std::string &name);
    
  private:
    
				     /**
				      * Abbreviate the somewhat
				      * lengthy name for the @p{Patch}
				      * class.
				      */
    typedef DataOutBase::Patch<2,2> Patch;

				     /**
				      * This is a list of patches that
				      * is created each time
				      * @p{build_patches} is
				      * called. These patches are used
				      * in the output routines of the
				      * base classes.
				      */
    std::vector<Patch> patches;

				     /**
				      * Name of the matrix to be
				      * written.
				      */
    std::string name;
    
				     /**
				      * Function by which the base
				      * class's functions get to know
				      * what patches they shall write
				      * to a file.
				      */
    virtual const std::vector<Patch> &
    get_patches () const;

				     /**
				      * Virtual function through which
				      * the names of data sets are
				      * obtained by the output
				      * functions of the base class.
				      */
    virtual std::vector<std::string> get_dataset_names () const;

				     /**
				      * Return the element with given
				      * indices of a sparse matrix.
				      */
    template <typename number>
    static double get_element (const SparseMatrix<number> &matrix,
			       const unsigned int          i,
			       const unsigned int          j);
    
				     /**
				      * Return the element with given
				      * indices of a block sparse
				      * matrix.
				      */
    template <typename number>
    static double get_element (const BlockSparseMatrix<number> &matrix,
			       const unsigned int               i,
			       const unsigned int               j);
    
				     /**
				      * Return the element with given
				      * indices from any matrix type
				      * for which no specialization of
				      * this function was declared
				      * above. This will call
				      * @p{operator()} on the matrix.
				      */
    template <class Matrix>
    static double get_element (const Matrix       &matrix,
			       const unsigned int  i,
			       const unsigned int  j);
    
};




/* ---------------------- Template and inline functions ------------- */


template <typename number>
inline
double
MatrixOut::get_element (const SparseMatrix<number> &matrix,
			const unsigned int          i,
			const unsigned int          j)
{
  return matrix.el(i,j);
};




template <typename number>
inline
double
MatrixOut::get_element (const BlockSparseMatrix<number> &matrix,
			const unsigned int               i,
			const unsigned int               j)
{
  return matrix.el(i,j);
};




template <class Matrix>
inline
double
MatrixOut::get_element (const Matrix       &matrix,
			const unsigned int  i,
			const unsigned int  j)
{
  return matrix(i,j);
};



template <class Matrix>
void
MatrixOut::build_patches (const Matrix      &matrix,
			  const std::string &name)
{
				   // first clear old data and set it
				   // to virgin state
  patches.clear ();
  patches.resize ((matrix.n()-1) * (matrix.m()-1));

				   // now build the patches
  unsigned int index=0;
  for (unsigned int i=0; i<matrix.m()-1; ++i)
    for (unsigned int j=0; j<matrix.n()-1; ++j, ++index)
      {
					 // within each patch, order
					 // the points in such a way
					 // that if some graphical
					 // output program (such as
					 // gnuplot) plots the
					 // quadrilaterals as two
					 // triangles, then the
					 // diagonal of the
					 // quadrilateral which cuts
					 // it into the two printed
					 // triangles is parallel to
					 // the diagonal of the
					 // matrix, rather than
					 // perpendicular to it. this
					 // has the advantage that,
					 // for example, the unit
					 // matrix is plotted as a
					 // straight rim, rather than
					 // as a series of bumps and
					 // valleys along the diagonal
	patches[index].vertices[0](0) = j;
	patches[index].vertices[0](1) = static_cast<signed int>(-i);
	patches[index].vertices[1](0) = j;
	patches[index].vertices[1](1) = static_cast<signed int>(-i-1);
	patches[index].vertices[2](0) = j+1;
	patches[index].vertices[2](1) = static_cast<signed int>(-i-1);
	patches[index].vertices[3](0) = j+1;
	patches[index].vertices[3](1) = static_cast<signed int>(-i);

	patches[index].n_subdivisions = 1;

	patches[index].data.reinit (1,4);
	patches[index].data(0,0) = get_element(matrix, i, j);
	patches[index].data(0,1) = get_element(matrix, i, j+1);
	patches[index].data(0,2) = get_element(matrix, i+1, j);
	patches[index].data(0,3) = get_element(matrix, i+1, j+1);
      };

				   // finally set the name
  this->name = name;
};



/*----------------------------   matrix_out.h     ---------------------------*/

#endif
/*----------------------------   matrix_out.h     ---------------------------*/


