// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__matrix_out_h
#define __deal2__matrix_out_h

#include <deal.II/base/config.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

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
 * @code
 *    FullMatrix<double> M;
 *    ...                // fill matrix M with some values
 *
 *                       // now write out M:
 *    MatrixOut matrix_out;
 *    std::ofstream out ("M.gnuplot");
 *    matrix_out.build_patches (M, "M");
 *    matrix_out.write_gnuplot (out);
 * @endcode
 * Of course, you can as well choose a different graphical output
 * format. Also, this class supports any matrix, not only of type
 * FullMatrix, as long as it satisfies a number of requirements,
 * stated with the member functions of this class.
 *
 * The generation of patches through the build_patches() function
 * can be modified by giving it an object holding certain flags. See
 * the documentation of the members of the Options class for a
 * description of these flags.
 *
 *
 * <h3>Internals</h3>
 *
 * To avoid a compiler error in Sun's Forte compiler, we derive
 * privately from DataOutBase. Since the base class
 * DataOutInterface does so as well, this does no harm, but
 * calms the compiler which is suspecting an access control conflict
 * otherwise. Testcase here:
 * @code
 *    template <typename T> class V {};
 *
 *    struct B1 {
 *        template <int dim> struct X {
 *      int i[dim];
 *        };
 *    };
 *
 *    struct B2 : private B1 {};
 *
 *    struct D : public B2, private B1 {
 *        ~D () {};
 *        typedef B1::X<2> X;
 *        V<X> x;
 *    };
 *
 *    D d;
 * @endcode
 *
 * @ingroup output
 * @author Wolfgang Bangerth, 2001
 */
class MatrixOut : public DataOutInterface<2,2>
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Class holding various
   * variables which are used to
   * modify the output of the
   * MatrixOut class.
   */
  struct Options
  {
    /**
     * If @p true, only show the
     * absolute values of the
     * matrix entries, rather
     * than their true values
     * including the
     * sign. Default value is
     * @p false.
     */
    bool         show_absolute_values;

    /**
     * If larger than one, do not
     * show each element of the
     * matrix, but rather an
     * average over a number of
     * entries. The number of
     * output patches is
     * accordingly smaller. This
     * flag determines how large
     * each shown block shall be
     * (in rows/columns). For
     * example, if it is two,
     * then always four entries
     * are collated into one.
     *
     * Default value is one.
     */
    unsigned int block_size;

    /**
     * If true, plot
     * discontinuous patches, one
     * for each entry.
     */
    bool discontinuous;

    /**
     * Default constructor. Set
     * all elements of this
     * structure to their default
     * values.
     */
    Options (const bool         show_absolute_values = false,
             const unsigned int block_size           = 1,
             const bool         discontinuous        = false);
  };

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
   * You may give a structure
   * holding various options. See
   * the description of the fields
   * of this structure for more
   * information.
   *
   * Note that this function
   * requires that we can extract
   * elements of the matrix, which
   * is done using the
   * get_element() function
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
   * <tt>m()</tt> and <tt>n()</tt>, which
   * return the number of rows and
   * columns, respectively.
   */
  template <class Matrix>
  void build_patches (const Matrix      &matrix,
                      const std::string &name,
                      const Options      options = Options(false, 1, false));

private:

  /**
   * Abbreviate the somewhat
   * lengthy name for the dealii::DataOutBase::Patch
   * class.
   *
   * Note that we have to indicate the
   * global scope in front of DataOutBase,
   * since otherwise the C++ rules specify
   * that this here indicates the
   * DataOutBase base class of this
   * class. Since that is a private base
   * class, we cannot access its members,
   * and so access to the local Patch type
   * would be forbidden.
   */
  typedef dealii::DataOutBase::Patch<2,2> Patch;

  /**
   * This is a list of patches that
   * is created each time
   * build_patches() is
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
                             const size_type             i,
                             const size_type             j);

  /**
   * Return the element with given
   * indices of a block sparse
   * matrix.
   */
  template <typename number>
  static double get_element (const BlockSparseMatrix<number> &matrix,
                             const size_type                  i,
                             const size_type                  j);

  /**
   * Return the element with given
   * indices from any matrix type
   * for which no specialization of
   * this function was declared
   * above. This will call
   * <tt>operator()</tt> on the matrix.
   */
  template <class Matrix>
  static double get_element (const Matrix       &matrix,
                             const size_type     i,
                             const size_type     j);

  /**
   * Get the value of the matrix at
   * gridpoint <tt>(i,j)</tt>. Depending
   * on the given flags, this can
   * mean different things, for
   * example if only absolute
   * values shall be shown then the
   * absolute value of the matrix
   * entry is taken. If the block
   * size is larger than one, then
   * an average of several matrix
   * entries is taken.
   */
  template <class Matrix>
  static double get_gridpoint_value (const Matrix       &matrix,
                                     const size_type     i,
                                     const size_type     j,
                                     const Options      &options);
};


/* ---------------------- Template and inline functions ------------- */


template <typename number>
inline
double
MatrixOut::get_element (const SparseMatrix<number> &matrix,
                        const size_type             i,
                        const size_type             j)
{
  return matrix.el(i,j);
}




template <typename number>
inline
double
MatrixOut::get_element (const BlockSparseMatrix<number> &matrix,
                        const size_type                  i,
                        const size_type                  j)
{
  return matrix.el(i,j);
}




template <class Matrix>
inline
double
MatrixOut::get_element (const Matrix   &matrix,
                        const size_type i,
                        const size_type j)
{
  return matrix(i,j);
}



template <class Matrix>
inline
double
MatrixOut::get_gridpoint_value (const Matrix   &matrix,
                                const size_type i,
                                const size_type j,
                                const Options  &options)
{
  // special case if block size is
  // one since we then don't need all
  // that loop overhead
  if (options.block_size == 1)
    {
      if (options.show_absolute_values == true)
        return std::fabs(get_element (matrix, i, j));
      else
        return get_element (matrix, i, j);
    }

  // if blocksize greater than one,
  // then compute average of elements
  double average = 0;
  size_type n_elements = 0;
  for (size_type row=i*options.block_size;
       row < std::min(size_type(matrix.m()),
                      size_type((i+1)*options.block_size)); ++row)
    for (size_type col=j*options.block_size;
         col < std::min(size_type(matrix.m()),
                        size_type((j+1)*options.block_size)); ++col, ++n_elements)
      if (options.show_absolute_values == true)
        average += std::fabs(get_element (matrix, row, col));
      else
        average += get_element (matrix, row, col);
  average /= n_elements;
  return average;
}



template <class Matrix>
void
MatrixOut::build_patches (const Matrix      &matrix,
                          const std::string &name,
                          const Options      options)
{
  size_type
  gridpoints_x = (matrix.n() / options.block_size
                  +
                  (matrix.n() % options.block_size != 0 ? 1 : 0)),
                 gridpoints_y = (matrix.m() / options.block_size
                                 +
                                 (matrix.m() % options.block_size != 0 ? 1 : 0));

  // If continuous, the number of
  // plotted patches is matrix size-1
  if (!options.discontinuous)
    {
      --gridpoints_x;
      --gridpoints_y;
    }

  // first clear old data and set it
  // to virgin state
  patches.clear ();
  patches.resize ((gridpoints_x) * (gridpoints_y));

  // now build the patches
  size_type index=0;
  for (size_type i=0; i<gridpoints_y; ++i)
    for (size_type j=0; j<gridpoints_x; ++j, ++index)
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
        patches[index].vertices[2](1) = static_cast<signed int>(-i);
        patches[index].vertices[3](0) = j+1;
        patches[index].vertices[3](1) = static_cast<signed int>(-i-1);
        // next scale all the patch
        // coordinates by the block
        // size, to get original
        // coordinates
        for (unsigned int v=0; v<4; ++v)
          patches[index].vertices[v] *= options.block_size;

        patches[index].n_subdivisions = 1;

        patches[index].data.reinit (1,4);
        if (options.discontinuous)
          {
            patches[index].data(0,0) = get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0,1) = get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0,2) = get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0,3) = get_gridpoint_value(matrix, i, j, options);
          }
        else
          {
            patches[index].data(0,0) = get_gridpoint_value(matrix, i,   j,   options);
            patches[index].data(0,1) = get_gridpoint_value(matrix, i+1, j,   options);
            patches[index].data(0,2) = get_gridpoint_value(matrix, i,   j+1, options);
            patches[index].data(0,3) = get_gridpoint_value(matrix, i+1, j+1, options);
          }
      };

  // finally set the name
  this->name = name;
}



/*----------------------------   matrix_out.h     ---------------------------*/

DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------   matrix_out.h     ---------------------------*/
