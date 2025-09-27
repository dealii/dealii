// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_matrix_out_h
#  define dealii_matrix_out_h

#  include <deal.II/base/config.h>

#  include <deal.II/base/data_out_base.h>

#  include <deal.II/lac/block_sparse_matrix.h>
#  include <deal.II/lac/sparse_matrix.h>

#  ifdef DEAL_II_WITH_TRILINOS
#    include <deal.II/lac/trilinos_block_sparse_matrix.h>
#    include <deal.II/lac/trilinos_sparse_matrix.h>
#  endif


DEAL_II_NAMESPACE_OPEN

/**
 * Output a matrix in graphical form using the generic format independent
 * output routines of the base class. The matrix is converted into a list of
 * patches on a 2d domain where the height is given by the elements of the
 * matrix. The functions of the base class can then write this "mountain
 * representation" of the matrix in a variety of graphical output formats. The
 * coordinates of the matrix output are that the columns run with increasing
 * x-axis, as usual, starting from zero, while the rows run into the negative
 * y-axis, also starting from zero. Note that due to some internal
 * restrictions, this class can only output one matrix at a time, i.e. it can
 * not take advantage of the multiple dataset capabilities of the base class.
 *
 * A typical usage of this class would be as follows:
 * @code
 *   FullMatrix<double> M;
 *   // fill matrix M with some values
 *   ...
 *
 *   // now write out M:
 *   MatrixOut matrix_out;
 *   std::ofstream out ("M.vtu");
 *   matrix_out.build_patches (M, "M");
 *   matrix_out.write_vtu (out);
 * @endcode
 * Of course, you can as well choose a different graphical output format:
 * all of the formats implemented by the DataOutInterface functions are also
 * supported by this class.
 *
 * This class supports any matrix, not only of type FullMatrix, as long
 * as it satisfies a number of requirements, stated with the member functions
 * of this class.
 *
 * The generation of patches through the build_patches() function can be
 * modified by giving it an object holding certain flags. See the
 * documentation of the members of the Options class for a description of
 * these flags.
 *
 *
 * @ingroup output
 */
class MatrixOut : public DataOutInterface<2, 2>
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Class holding various variables which are used to modify the output of
   * the MatrixOut class.
   */
  struct Options
  {
    /**
     * If @p true, only show the absolute values of the matrix entries, rather
     * than their true values including the sign. This is useful if you have
     * matrix entries that span a large range and want to display their
     * magnitude using a logarithmic scale.
     *
     * The default value is @p false.
     */
    bool show_absolute_values;

    /**
     * If larger than one, do not show each element of the matrix, but rather
     * an average over a number of entries. The number of output patches is
     * accordingly smaller. This flag determines how large each shown block
     * shall be (in rows/columns). For example, if it is two, then always four
     * entries are collated into one.
     *
     * The default value is one.
     */
    unsigned int block_size;

    /**
     * By default, this class shows each matrix entry as one value,
     * and then produces a bilinear display of all of these
     * values. This results in a continuous plot on a mesh with
     * $N\times N$ vertices for a matrix of size $N$. If, on the other
     * hand, the current variable is set to `true`, then this class
     * instead creates a plot in which each matrix entry corresponds
     * to a single patch on which the value of the matrix entry is
     * shown as a constant function. This creates a discontinuous plot
     * with $N\times N$ cells on a mesh of size $(N+1)\times(N+1)$.
     *
     * The default value is false.
     */
    bool discontinuous;

    /**
     * If this flag is set to `false`, the MatrixOut class creates
     * a plot in which each matrix
     * entry is actually shown, even if it is zero. For large
     * matrices, this results in very large output files, or indeed
     * exhausts the available memory. On the other hand, if the
     * current flags is set to `true`, then the class only outputs
     * patches whenever there are non-zero matrix entries to be
     * shown. For sparse matrices, this leads to an output size that
     * is proportional to the number of nonzero entries, rather than
     * proportional to $N^2$.
     *
     * The default is `true`.
     *
     * @note Internally, the current implementation continues to loop
     *   over all matrix entries, whether they are zero or not. As a
     *   consequence, the *run time* of outputting sparse matrices
     *   continues to be proportional to $N^2$, rather than proportional
     *   to the number of nonzero entries. For large matrices, this means
     *   that you still have to have patience -- but at least it is possible
     *   to output information about matrices of size $10,000\times 10,000$
     *   or $50,000\times 50,000$ with a few million nonzero entries.
     */
    bool create_sparse_plot;

    /**
     * Default constructor. Set all elements of this structure to their
     * default values.
     */
    Options(const bool         show_absolute_values = false,
            const unsigned int block_size           = 1,
            const bool         discontinuous        = false,
            const bool         create_sparse_plot   = true);
  };

  /**
   * Destructor. Declared in order to make it virtual.
   */
  virtual ~MatrixOut() override = default;

  /**
   * Generate a list of patches from the given matrix and use the given string
   * as the name of the data set upon writing to a file. Once patches have
   * been built, you can use the functions of the base class to write the data
   * into a files, using one of the supported output formats.
   *
   * The last argument provides customization choices for how output
   * is to be produced. See the description of the fields of the
   * Options structure for more information.
   *
   * Note that this function requires that we can extract elements of the
   * matrix, which is done using the get_element() function declared in an
   * internal namespace. By adding specializations, you can extend this class
   * to other matrix classes which are not presently supported. Furthermore,
   * we need to be able to extract the size of the matrix, for which we assume
   * that the matrix type offers member functions <tt>m()</tt> and
   * <tt>n()</tt>, which return the number of rows and columns, respectively.
   */
  template <class Matrix>
  void
  build_patches(const Matrix      &matrix,
                const std::string &name,
                const Options      options = Options());

private:
  /**
   * Abbreviate the somewhat lengthy name for the dealii::DataOutBase::Patch
   * class.
   */
  using Patch = DataOutBase::Patch<2, 2>;

  /**
   * This is a list of patches that is created each time build_patches() is
   * called. These patches are used in the output routines of the base
   * classes.
   */
  std::vector<Patch> patches;

  /**
   * Name of the matrix to be written.
   */
  std::string name;

  /**
   * %Function by which the base class's functions get to know what patches
   * they shall write to a file.
   */
  virtual const std::vector<Patch> &
  get_patches() const override;

  /**
   * Virtual function through which the names of data sets are obtained by the
   * output functions of the base class.
   */
  virtual std::vector<std::string>
  get_dataset_names() const override;

  /**
   * Get the value of the matrix at gridpoint <tt>(i,j)</tt>. Depending on the
   * given flags, this can mean different things, for example if only absolute
   * values shall be shown then the absolute value of the matrix entry is
   * taken. If the block size is larger than one, then an average of several
   * matrix entries is taken.
   */
  template <class Matrix>
  static double
  get_gridpoint_value(const Matrix   &matrix,
                      const size_type i,
                      const size_type j,
                      const Options  &options);
};


/* ---------------------- Template and inline functions ------------- */


namespace internal
{
  namespace MatrixOutImplementation
  {
    /**
     * Return the element with given indices of a sparse matrix.
     */
    template <typename number>
    double
    get_element(const dealii::SparseMatrix<number> &matrix,
                const types::global_dof_index       i,
                const types::global_dof_index       j)
    {
      return matrix.el(i, j);
    }



    /**
     * Return the element with given indices of a block sparse matrix.
     */
    template <typename number>
    double
    get_element(const dealii::BlockSparseMatrix<number> &matrix,
                const types::global_dof_index            i,
                const types::global_dof_index            j)
    {
      return matrix.el(i, j);
    }


#  ifdef DEAL_II_WITH_TRILINOS
    /**
     * Return the element with given indices of a Trilinos sparse matrix.
     */
    inline double
    get_element(const TrilinosWrappers::SparseMatrix &matrix,
                const types::global_dof_index         i,
                const types::global_dof_index         j)
    {
      return matrix.el(i, j);
    }



    /**
     * Return the element with given indices of a Trilinos block sparse
     * matrix.
     */
    inline double
    get_element(const TrilinosWrappers::BlockSparseMatrix &matrix,
                const types::global_dof_index              i,
                const types::global_dof_index              j)
    {
      return matrix.el(i, j);
    }
#  endif


#  ifdef DEAL_II_WITH_PETSC
    // no need to do anything: PETSc matrix objects do not distinguish
    // between operator() and el(i,j), so we can safely access elements
    // through the generic function below
#  endif


    /**
     * Return the element with given indices from any matrix type for which
     * no specialization of this function was declared above. This will call
     * <tt>operator()</tt> on the matrix.
     */
    template <class Matrix>
    double
    get_element(const Matrix                 &matrix,
                const types::global_dof_index i,
                const types::global_dof_index j)
    {
      return matrix(i, j);
    }
  } // namespace MatrixOutImplementation
} // namespace internal



template <class Matrix>
inline double
MatrixOut::get_gridpoint_value(const Matrix   &matrix,
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
        return std::fabs(
          internal::MatrixOutImplementation::get_element(matrix, i, j));
      else
        return internal::MatrixOutImplementation::get_element(matrix, i, j);
    }

  // if blocksize greater than one,
  // then compute average of elements
  double    average    = 0;
  size_type n_elements = 0;
  for (size_type row = i * options.block_size;
       row <
       std::min(size_type(matrix.m()), size_type((i + 1) * options.block_size));
       ++row)
    for (size_type col = j * options.block_size;
         col < std::min(size_type(matrix.m()),
                        size_type((j + 1) * options.block_size));
         ++col, ++n_elements)
      if (options.show_absolute_values == true)
        average += std::fabs(
          internal::MatrixOutImplementation::get_element(matrix, row, col));
      else
        average +=
          internal::MatrixOutImplementation::get_element(matrix, row, col);
  average /= n_elements;
  return average;
}



template <class Matrix>
void
MatrixOut::build_patches(const Matrix      &matrix,
                         const std::string &name,
                         const Options      options)
{
  size_type n_patches_x = (matrix.n() / options.block_size +
                           (matrix.n() % options.block_size != 0 ? 1 : 0)),
            n_patches_y = (matrix.m() / options.block_size +
                           (matrix.m() % options.block_size != 0 ? 1 : 0));

  // If continuous, the number of
  // plotted patches is matrix size-1
  if (!options.discontinuous)
    {
      --n_patches_x;
      --n_patches_y;
    }

  const size_type n_patches =
    (options.create_sparse_plot ?
       [&]() {
         size_type count = 0;
         for (size_type i = 0; i < n_patches_y; ++i)
           {
             for (size_type j = 0; j < n_patches_x; ++j)
               // Use the same logic as below to determine whether we
               // need to output a patch, and count if we do:
               if ((((options.discontinuous == true) &&
                     (get_gridpoint_value(matrix, i, j, options) != 0)) ||
                    ((options.discontinuous == false) &&
                     ((get_gridpoint_value(matrix, i, j, options) != 0) ||
                      (get_gridpoint_value(matrix, i + 1, j, options) != 0) ||
                      (get_gridpoint_value(matrix, i, j + 1, options) != 0) ||
                      (get_gridpoint_value(matrix, i + 1, j + 1, options) !=
                       0)))))
                 ++count;
           }
         return count;
       }() :
       n_patches_x * n_patches_y);

  // first clear old data and re-set the object to a correctly sized state:
  patches.clear();
  try
    {
      patches.resize(n_patches);
    }
  catch (const std::bad_alloc &)
    {
      AssertThrow(false,
                  ExcMessage("You are trying to create a graphical "
                             "representation of a matrix that would "
                             "requiring outputting " +
                             (options.create_sparse_plot ?
                                std::to_string(n_patches) :
                                std::to_string(n_patches_x) + "x" +
                                  std::to_string(n_patches_y)) +
                             " patches. There is not enough memory to output " +
                             "this many patches."));
    }

  // now build the patches
  size_type index = 0;
  for (size_type i = 0; i < n_patches_y; ++i)
    for (size_type j = 0; j < n_patches_x; ++j)
      {
        // If we are creating a sparse plot, check whether this patch
        // would have any nonzero values. If not, we can skip the
        // patch:
        if (options.create_sparse_plot &&
            (((options.discontinuous == true) &&
              (get_gridpoint_value(matrix, i, j, options) == 0)) ||
             ((options.discontinuous == false) &&
              (get_gridpoint_value(matrix, i, j, options) == 0) &&
              (get_gridpoint_value(matrix, i + 1, j, options) == 0) &&
              (get_gridpoint_value(matrix, i, j + 1, options) == 0) &&
              (get_gridpoint_value(matrix, i + 1, j + 1, options) == 0))))
          continue;

        patches[index].n_subdivisions = 1;
        patches[index].reference_cell = ReferenceCells::Quadrilateral;

        // within each patch, order the points in such a way that if some
        // graphical output program (such as gnuplot) plots the quadrilaterals
        // as two triangles, then the diagonal of the quadrilateral which cuts
        // it into the two printed triangles is parallel to the diagonal of the
        // matrix, rather than perpendicular to it. this has the advantage that,
        // for example, the unit matrix is plotted as a straight ridge, rather
        // than as a series of bumps and valleys along the diagonal
        patches[index].vertices[0][0] = j;
        patches[index].vertices[0][1] = -static_cast<signed int>(i);
        patches[index].vertices[1][0] = j;
        patches[index].vertices[1][1] = -static_cast<signed int>(i + 1);
        patches[index].vertices[2][0] = j + 1;
        patches[index].vertices[2][1] = -static_cast<signed int>(i);
        patches[index].vertices[3][0] = j + 1;
        patches[index].vertices[3][1] = -static_cast<signed int>(i + 1);
        // next scale all the patch
        // coordinates by the block
        // size, to get original
        // coordinates
        for (auto &vertex : patches[index].vertices)
          vertex *= options.block_size;

        patches[index].n_subdivisions = 1;

        patches[index].data.reinit(1, 4);
        if (options.discontinuous)
          {
            patches[index].data(0, 0) =
              get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0, 1) =
              get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0, 2) =
              get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0, 3) =
              get_gridpoint_value(matrix, i, j, options);
          }
        else
          {
            patches[index].data(0, 0) =
              get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0, 1) =
              get_gridpoint_value(matrix, i + 1, j, options);
            patches[index].data(0, 2) =
              get_gridpoint_value(matrix, i, j + 1, options);
            patches[index].data(0, 3) =
              get_gridpoint_value(matrix, i + 1, j + 1, options);
          }

        ++index;
      }

  // finally set the name
  this->name = name;
}



/*----------------------------   matrix_out.h     ---------------------------*/

DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------   matrix_out.h     ---------------------------*/
