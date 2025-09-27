// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tensor_product_matrix_h
#define dealii_tensor_product_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/floating_point_comparator.h>
#include <deal.II/base/mutex.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/matrix_free/tensor_product_kernels.h>

#include <bitset>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename>
class Vector;
template <typename>
class FullMatrix;
#endif

/**
 * This is a special matrix class defined as the tensor product (or Kronecker
 * product) of 1d matrices of the type
 * @f{align*}{
 * L &= A_1 \otimes M_0 + M_1 \otimes A_0
 * @f}
 * in 2d and
 * @f{align*}{
 * L &= A_2 \otimes M_1 \otimes M_0 + M_2 \otimes A_1 \otimes M_0 + M_2 \otimes
 * M_1 \otimes A_0
 * @f}
 * in 3d. The typical application setting is a discretization of the Laplacian
 * $L$ on a Cartesian (axis-aligned) geometry, where it can be exactly
 * represented by the Kronecker or tensor product of a 1d @ref GlossMassMatrix "mass matrix" $M$ and
 * a 1d Laplace matrix $A$ in each tensor direction (due to symmetry $M$ and $A$
 * are the same in each dimension). The dimension of the resulting class is the
 * product of the one-dimensional matrices.
 *
 * This class implements two basic operations, namely the usual multiplication
 * by a vector and the inverse. For both operations, fast tensorial techniques
 * can be applied that implement the operator evaluation in
 * $\text{size}(M)^{d+1}$ arithmetic operations, considerably less than
 * $\text{size}(M)^{2d}$ for the naive forward transformation and
 * $\text{size}(M)^{3d}$ for setting up the inverse of $L$.
 *
 * Interestingly, the exact inverse of the matrix $L$ can be found through
 * tensor products due to 1964's work by Lynch et al. @cite Lynch1964,
 * @f{align*}{
 * L^{-1} &= S_1 \otimes S_0 (\Lambda_1 \otimes I + I \otimes \Lambda_0)^{-1}
 * S_1^\mathrm T \otimes S_0^\mathrm T,
 * @f}
 * where $S_d$ is the matrix of eigenvectors to the generalized eigenvalue
 * problem in the given tensor direction $d$:
 * @f{align*}{
 * A_d s  &= \lambda M_d s, d = 0, \quad \ldots,\mathrm{dim},
 * @f}
 * and $\Lambda_d$ is the diagonal matrix representing the generalized
 * eigenvalues $\lambda$. Note that the vectors $s$ are such that they
 * simultaneously diagonalize $A_d$ and $M_d$, i.e. $S_d^{\mathrm T} A_d S_d =
 * \Lambda_d$ and $S_d^{\mathrm T} M_d S_d = I$. This method of matrix inversion
 * is called fast diagonalization method.
 *
 * This class requires LAPACK support.
 *
 * @note This class allows for two modes of usage. The first is a use case
 * with run time constants for the matrix dimensions that is achieved by
 * setting the optional template parameter <tt>n_rows_1d</tt> to -1. The second
 * mode of usage that is faster allows to set the template parameter as a
 * compile time constant, giving significantly faster code in particular for
 * small sizes of the matrix.
 *
 * @note This class can work with scalar types (float, double) and
 * VectorizedArray types.
 *
 * @tparam dim Dimension of the problem. Currently, 1D, 2d, and 3d codes are
 * implemented.
 *
 * @tparam Number Arithmetic type of the underlying array elements. Note that the
 * underlying LAPACK implementation supports only float and double numbers, so
 * only these two types are currently supported by the generic class.
 * Nevertheless, a template specialization for the vectorized types
 * VectorizedArray<float> and VectorizedArray<double> exists. This is necessary
 * to perform LAPACK calculations for each vectorization lane, i.e. for the
 * supported float and double numbers.
 *
 * @tparam n_rows_1d Compile-time number of rows of 1d matrices (only
 * valid if the number of rows and columns coincide for each
 * dimension). By default at -1, which means that the number of rows
 * is determined at run-time by means of the matrices passed to the
 * reinit() function.
 */
template <int dim, typename Number, int n_rows_1d = -1>
class TensorProductMatrixSymmetricSum
{
public:
  /**
   * Type of matrix entries. This alias is analogous to <tt>value_type</tt>
   * in the standard library containers.
   */
  using value_type = Number;

  /**
   * The static number of rows of the 1d matrices. For more details,
   * see the description of the template parameter <tt>n_rows_1d</tt>.
   */
  static constexpr int n_rows_1d_static = n_rows_1d;

  /**
   * Default constructor.
   */
  TensorProductMatrixSymmetricSum() = default;

  /**
   * Constructor that is equivalent to the empty constructor and
   * immediately calling reinit(mass_matrix, derivative_matrix).
   */
  template <typename T>
  TensorProductMatrixSymmetricSum(const T &mass_matrix,
                                  const T &derivative_matrix);

  /**
   * Initializes the tensor product matrix by copying the arrays of 1d mass
   * matrices @p mass_matrix and 1d derivative matrices @p derivative_matrix into its
   * base class counterparts, respectively, and by assembling the regarding
   * generalized eigenvalues and eigenvectors in eigenvalues
   * and eigenvectors, respectively.
   * Note that the current implementation requires each $M_{d}$ to be symmetric
   * and positive definite and every $A_{d}$ to be symmetric and invertible but
   * not necessarily positive definite. Columns and rows filled with zero are
   * ignored.
   *
   * @warning This class accepts the following types:
   * "std::array<Table<2, Number>, dim>", "std::array<FullMatrix<Number>, dim>",
   * and "Table<2, Number>". In the latter case, we consider the same 1d
   * @ref GlossMassMatrix "mass matrix" @p mass_matrix and the same 1d derivative matrix
   * @p derivative_matrix for each tensor direction.
   */
  template <typename T>
  void
  reinit(const T &mass_matrix, const T &derivative_matrix);

  /**
   * Return the number of rows of the tensor product matrix
   * resulting from the Kronecker product of 1d matrices, which is described
   * in the main documentation of TensorProductMatrixSymmetricSum.
   */
  unsigned int
  m() const;

  /**
   * Return the number of columns of the tensor product matrix
   * resulting from the Kronecker product of 1d matrices, which is described
   * in the main documentation of TensorProductMatrixSymmetricSum.
   */
  unsigned int
  n() const;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of TensorProductMatrixSymmetricSum.
   * This function is operating on ArrayView to allow checks of
   * array bounds with respect to @p dst and @p src.
   *
   * @warning This function works on an internal temporal array, leading to
   * increased memory consumption if many instances of this class are created,
   * e.g., a different object on every cell with different underlying
   * coefficients each. Furthermore, only one thread runs this function at a
   * time (ensured internally with a mutex). If these two limitations are an
   * issue, please consider the other version of this function.
   */
  void
  vmult(const ArrayView<Number> &dst, const ArrayView<const Number> &src) const;

  /**
   * Same as above but letting the user provide a user-owned temporary array,
   * resolving the two issues described above. This array is resized
   * internally to the needed size.
   */
  void
  vmult(const ArrayView<Number>       &dst,
        const ArrayView<const Number> &src,
        AlignedVector<Number>         &tmp) const;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of TensorProductMatrixSymmetricSum.
   * This function is operating on ArrayView to allow checks of
   * array bounds with respect to @p dst and @p src.
   */
  void
  apply_inverse(const ArrayView<Number>       &dst,
                const ArrayView<const Number> &src) const;

  /**
   * Return the memory consumption of the allocated memory in this class.
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * An array containing a @ref GlossMassMatrix "mass matrix" for each tensor direction.
   */
  std::array<Table<2, Number>, dim> mass_matrix;

  /**
   * An array containing a derivative matrix for each tensor direction.
   */
  std::array<Table<2, Number>, dim> derivative_matrix;

  /**
   * An array storing the generalized eigenvalues
   * for each tensor direction.
   */
  std::array<AlignedVector<Number>, dim> eigenvalues;

  /**
   * An array storing the generalized eigenvectors
   * for each tensor direction.
   */
  std::array<Table<2, Number>, dim> eigenvectors;

private:
  /**
   * An array for temporary data.
   */
  mutable AlignedVector<Number> tmp_array;

  /**
   * A mutex that guards access to the array @p tmp_array.
   */
  mutable Threads::Mutex mutex;
};



namespace internal
{
  namespace TensorProductMatrixSymmetricSum
  {
    template <typename Number>
    struct MatrixPairComparator
    {
      using VectorizedArrayTrait =
        dealii::internal::VectorizedArrayTrait<Number>;
      using ScalarNumber = typename VectorizedArrayTrait::value_type;
      static constexpr std::size_t width = VectorizedArrayTrait::width();

      using MatrixPairType =
        std::pair<std::bitset<width>,
                  std::pair<Table<2, Number>, Table<2, Number>>>;

      MatrixPairComparator()
        : eps(std::sqrt(std::numeric_limits<ScalarNumber>::epsilon()))
      {}

      bool
      operator()(const MatrixPairType &left, const MatrixPairType &right) const
      {
        const auto &M_0 = left.second.first;
        const auto &K_0 = left.second.second;
        const auto &M_1 = right.second.first;
        const auto &K_1 = right.second.second;

        std::bitset<width> mask;

        for (unsigned int v = 0; v < width; ++v)
          mask[v] = left.first[v] && right.first[v];

        const FloatingPointComparator<Number> comparator(
          eps, false /*use relative tolerance*/, mask);

        if (comparator(M_0, M_1))
          return true;
        else if (comparator(M_1, M_0))
          return false;
        else if (comparator(K_0, K_1))
          return true;
        else
          return false;
      }

    private:
      const ScalarNumber eps;
    };
  } // namespace TensorProductMatrixSymmetricSum
} // namespace internal



/**
 * A class similar to TensorProductMatrixSymmetricSum.
 *
 * The class TensorProductMatrixSymmetricSum stores a
 * 1d @ref GlossMassMatrix "mass matrix", 1d @ref GlossStiffnessMatrix "stiffness matrix", eigenvalues and eigenvectors
 * for each direction. If one uses one TensorProductMatrixSymmetricSum
 * instance for, e.g., each cell, these quantities are stored
 * for each cell. There is no possibility to reuse quantities between
 * TensorProductMatrixSymmetricSum instances even if the values of the
 * internal data structures might be the same. This class targets the case
 * of many TensorProductMatrixSymmetricSum instances, where some of them might
 * possibly share the underlying 1d matrices and hence re-use the same data.
 *
 * This class is flexible and allows to interpret the parameter
 * @p index arbitrarily. In the case of an element-centric patch
 * smoother, the index might correspond to the cell index and,
 * in the case of a vertex patch smoother, the index might
 * correspond to the vertex index defining the vertex patch. If
 * @p n_rows_1d is set to -1, the sizes of the mass matrices and
 * the stiffness matrices can differ between cells, which might be
 * useful if the class is used in the context of hp-refinement to
 * construct a smoother.
 */
template <int dim, typename Number, int n_rows_1d = -1>
class TensorProductMatrixSymmetricSumCollection
{
  using MatrixPairType = std::pair<Table<2, Number>, Table<2, Number>>;

  using MatrixPairTypeWithMask = std::pair<
    std::bitset<dealii::internal::VectorizedArrayTrait<Number>::width()>,
    MatrixPairType>;

public:
  /**
   * Struct to configure TensorProductMatrixSymmetricSumCollection.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData(const bool compress_matrices           = true,
                   const bool precompute_inverse_diagonal = true);

    /**
     * Try to compress internal matrices. Default: true.
     */
    bool compress_matrices;

    /**
     * Precompute inverse diagonal.
     */
    bool precompute_inverse_diagonal;
  };

  /**
   * Constructor.
   */
  TensorProductMatrixSymmetricSumCollection(
    const AdditionalData &additional_data = AdditionalData());

  /**
   * Allocate memory. The parameter @p specifies the maximum value
   * of the index used in invert() and apply_inverse().
   */
  void
  reserve(const unsigned int size);

  /**
   * For a given @p index, attach the mass matrices @p Ms and
   * stiffness matrices @p Ks to the stored data, looking out for
   * compression possibilities.
   */
  template <typename T>
  void
  insert(const unsigned int index, const T &Ms, const T &Ks);

  /**
   * Finalize setup. This function computes, e.g., the
   * eigenvalues and the eigenvectors.
   */
  void
  finalize();

  /**
   * Apply the inverse matrix for a given @p index.
   */
  void
  apply_inverse(const unsigned int             index,
                const ArrayView<Number>       &dst_in,
                const ArrayView<const Number> &src_in) const;

  /**
   * Return the memory consumption of this class in bytes.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Return the number of 1d matrices of each type stored internally.
   * In the case that no compression could be performed, its value
   * is the parameter passed to the function reserve() times the
   * number of dimension. If compression could be performed, the
   * value returned is less. In the optimal case (uniform Cartesian
   * mesh), the value is one.
   */
  std::size_t
  storage_size() const;

private:
  /**
   * Try to compress matrices.
   */
  const bool compress_matrices;

  /**
   * Precompute inverse diagonal.
   */
  const bool precompute_inverse_diagonal;

  /**
   * Container used to collect 1d matrices if no compression is
   * requested. The memory is freed during finalize().
   */
  std::vector<MatrixPairType> mass_and_derivative_matrices;

  /**
   * Container used during setup to determine the unique 1d
   * matrices. The memory is freed during finalize().
   */
  std::map<
    MatrixPairTypeWithMask,
    unsigned int,
    internal::TensorProductMatrixSymmetricSum::MatrixPairComparator<Number>>
    cache;

  /**
   * Map from index to the storage position within mass_matrices,
   * derivative_matrices, eigenvectors, and
   * eigenvalues. If compression was not successful, this
   * vector is empty, since the vectors can be access directly
   * with the given index.
   */
  std::vector<unsigned int> indices;

  /**
   * Vector of 1d mass matrices.
   */
  AlignedVector<Number> mass_matrices;

  /**
   * Vector of 1d derivative matrices.
   */
  AlignedVector<Number> derivative_matrices;

  /**
   * Vector of eigenvectors.
   */
  AlignedVector<Number> eigenvectors;

  /**
   * Vector of eigenvalues.
   */
  AlignedVector<Number> eigenvalues;

  /**
   * Vector of inverted eigenvalues.
   */
  AlignedVector<Number> inverted_eigenvalues;

  /**
   * Pointer into mass_matrices, derivative_matrices, and eigenvalues.
   */
  std::vector<unsigned int> vector_ptr;

  /**
   * Pointer into mass_matrices, derivative_matrices, and eigenvalues.
   */
  std::vector<unsigned int> matrix_ptr;

  /**
   * Number of rows in 1 of each cell.
   */
  std::vector<unsigned int> vector_n_rows_1d;
};


/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN

namespace internal
{
  namespace TensorProductMatrixSymmetricSum
  {
    /**
     * Compute generalized eigenvalues and eigenvectors of the real
     * generalized symmetric eigenproblem $A v = \lambda M v$. Since we are
     * operating on plain pointers we require the size of the matrices
     * beforehand. Note that the data arrays for the eigenvalues and
     * eigenvectors have to be initialized to a proper size, too. (no check of
     * array bounds possible)
     */
    template <typename Number>
    void
    spectral_assembly(const Number      *mass_matrix,
                      const Number      *derivative_matrix,
                      const unsigned int n_rows,
                      const unsigned int n_cols,
                      Number            *eigenvalues,
                      Number            *eigenvectors)
    {
      Assert(n_rows == n_cols, ExcNotImplemented());

      std::vector<bool> constrained_dofs(n_rows, false);

      for (unsigned int i = 0; i < n_rows; ++i)
        {
          if (mass_matrix[i + i * n_rows] == 0.0)
            {
              Assert(derivative_matrix[i + i * n_rows] == 0.0,
                     ExcInternalError());

              for (unsigned int j = 0; j < n_rows; ++j)
                {
                  Assert(derivative_matrix[i + j * n_rows] == 0,
                         ExcInternalError());
                  Assert(derivative_matrix[j + i * n_rows] == 0,
                         ExcInternalError());
                }

              constrained_dofs[i] = true;
            }
        }

      const auto transpose_fill_nm = [&constrained_dofs](Number            *out,
                                                         const Number      *in,
                                                         const unsigned int n,
                                                         const unsigned int m) {
        for (unsigned int mm = 0, c = 0; mm < m; ++mm)
          for (unsigned int nn = 0; nn < n; ++nn, ++c)
            out[mm + nn * m] =
              (mm == nn && constrained_dofs[mm]) ? Number(1.0) : in[c];
      };

      std::vector<dealii::Vector<Number>> eigenvecs(n_rows);
      LAPACKFullMatrix<Number>            mass_copy(n_rows, n_cols);
      LAPACKFullMatrix<Number>            deriv_copy(n_rows, n_cols);

      transpose_fill_nm(&(mass_copy(0, 0)), mass_matrix, n_rows, n_cols);
      transpose_fill_nm(&(deriv_copy(0, 0)), derivative_matrix, n_rows, n_cols);

      deriv_copy.compute_generalized_eigenvalues_symmetric(mass_copy,
                                                           eigenvecs);
      AssertDimension(eigenvecs.size(), n_rows);
      for (unsigned int i = 0, c = 0; i < n_rows; ++i)
        for (unsigned int j = 0; j < n_cols; ++j, ++c)
          if (constrained_dofs[i] == false)
            eigenvectors[c] = eigenvecs[j][i];

      for (unsigned int i = 0; i < n_rows; ++i, ++eigenvalues)
        *eigenvalues = deriv_copy.eigenvalue(i).real();
    }



    template <std::size_t dim, typename Number>
    inline void
    setup(const std::array<Table<2, Number>, dim> &mass_matrix,
          const std::array<Table<2, Number>, dim> &derivative_matrix,
          std::array<Table<2, Number>, dim>       &eigenvectors,
          std::array<AlignedVector<Number>, dim>  &eigenvalues)
    {
      const unsigned int n_rows_1d = mass_matrix[0].n_cols();

      for (unsigned int dir = 0; dir < dim; ++dir)
        {
          AssertDimension(n_rows_1d, mass_matrix[dir].n_cols());
          AssertDimension(mass_matrix[dir].n_rows(), mass_matrix[dir].n_cols());
          AssertDimension(mass_matrix[dir].n_rows(),
                          derivative_matrix[dir].n_rows());
          AssertDimension(mass_matrix[dir].n_rows(),
                          derivative_matrix[dir].n_cols());

          eigenvectors[dir].reinit(mass_matrix[dir].n_cols(),
                                   mass_matrix[dir].n_rows());
          eigenvalues[dir].resize(mass_matrix[dir].n_cols());
          internal::TensorProductMatrixSymmetricSum::spectral_assembly<Number>(
            &(mass_matrix[dir](0, 0)),
            &(derivative_matrix[dir](0, 0)),
            mass_matrix[dir].n_rows(),
            mass_matrix[dir].n_cols(),
            eigenvalues[dir].begin(),
            &(eigenvectors[dir](0, 0)));
        }
    }



    template <std::size_t dim, typename Number, std::size_t n_lanes>
    inline void
    setup(
      const std::array<Table<2, VectorizedArray<Number, n_lanes>>, dim>
        &mass_matrix,
      const std::array<Table<2, VectorizedArray<Number, n_lanes>>, dim>
        &derivative_matrix,
      std::array<Table<2, VectorizedArray<Number, n_lanes>>, dim> &eigenvectors,
      std::array<AlignedVector<VectorizedArray<Number, n_lanes>>, dim>
        &eigenvalues)
    {
      const unsigned int     n_rows_1d = mass_matrix[0].n_cols();
      constexpr unsigned int macro_size =
        VectorizedArray<Number, n_lanes>::size();
      const std::size_t nm_flat_size_max = n_rows_1d * n_rows_1d * macro_size;
      const std::size_t n_flat_size_max  = n_rows_1d * macro_size;

      std::vector<Number> mass_matrix_flat;
      std::vector<Number> deriv_matrix_flat;
      std::vector<Number> eigenvalues_flat;
      std::vector<Number> eigenvectors_flat;
      mass_matrix_flat.resize(nm_flat_size_max);
      deriv_matrix_flat.resize(nm_flat_size_max);
      eigenvalues_flat.resize(n_flat_size_max);
      eigenvectors_flat.resize(nm_flat_size_max);
      std::array<unsigned int, macro_size> offsets_nm;
      std::array<unsigned int, macro_size> offsets_n;
      for (unsigned int dir = 0; dir < dim; ++dir)
        {
          AssertDimension(n_rows_1d, mass_matrix[dir].n_cols());
          AssertDimension(mass_matrix[dir].n_rows(), mass_matrix[dir].n_cols());
          AssertDimension(mass_matrix[dir].n_rows(),
                          derivative_matrix[dir].n_rows());
          AssertDimension(mass_matrix[dir].n_rows(),
                          derivative_matrix[dir].n_cols());

          const unsigned int n_rows = mass_matrix[dir].n_rows();
          const unsigned int n_cols = mass_matrix[dir].n_cols();
          const unsigned int nm     = n_rows * n_cols;
          for (unsigned int vv = 0; vv < macro_size; ++vv)
            offsets_nm[vv] = nm * vv;

          vectorized_transpose_and_store<Number, n_lanes>(
            false,
            nm,
            &(mass_matrix[dir](0, 0)),
            offsets_nm.data(),
            mass_matrix_flat.data());
          vectorized_transpose_and_store<Number, n_lanes>(
            false,
            nm,
            &(derivative_matrix[dir](0, 0)),
            offsets_nm.data(),
            deriv_matrix_flat.data());

          const Number *mass_cbegin    = mass_matrix_flat.data();
          const Number *deriv_cbegin   = deriv_matrix_flat.data();
          Number       *eigenvec_begin = eigenvectors_flat.data();
          Number       *eigenval_begin = eigenvalues_flat.data();
          for (unsigned int lane = 0; lane < macro_size; ++lane)
            internal::TensorProductMatrixSymmetricSum::spectral_assembly<
              Number>(mass_cbegin + nm * lane,
                      deriv_cbegin + nm * lane,
                      n_rows,
                      n_cols,
                      eigenval_begin + n_rows * lane,
                      eigenvec_begin + nm * lane);

          eigenvalues[dir].resize(n_rows);
          eigenvectors[dir].reinit(n_rows, n_cols);
          for (unsigned int vv = 0; vv < macro_size; ++vv)
            offsets_n[vv] = n_rows * vv;
          vectorized_load_and_transpose<Number, n_lanes>(
            n_rows,
            eigenvalues_flat.data(),
            offsets_n.data(),
            eigenvalues[dir].begin());
          vectorized_load_and_transpose<Number, n_lanes>(
            nm,
            eigenvectors_flat.data(),
            offsets_nm.data(),
            &(eigenvectors[dir](0, 0)));
        }
    }



    template <std::size_t dim, typename Number>
    inline std::array<Table<2, Number>, dim>
    convert(const std::array<Table<2, Number>, dim> &mass_matrix)
    {
      return mass_matrix;
    }



    template <std::size_t dim, typename Number>
    inline std::array<Table<2, Number>, dim>
    convert(const std::array<FullMatrix<Number>, dim> &mass_matrix)
    {
      std::array<Table<2, Number>, dim> mass_copy;

      std::transform(mass_matrix.cbegin(),
                     mass_matrix.cend(),
                     mass_copy.begin(),
                     [](const FullMatrix<Number> &m) -> Table<2, Number> {
                       return m;
                     });

      return mass_copy;
    }



    template <std::size_t dim, typename Number>
    inline std::array<Table<2, Number>, dim>
    convert(const Table<2, Number> &matrix)
    {
      std::array<Table<2, Number>, dim> matrices;

      std::fill(matrices.begin(), matrices.end(), matrix);

      return matrices;
    }



    template <int n_rows_1d_templated, std::size_t dim, typename Number>
    void
    vmult(Number                                *dst,
          const Number                          *src,
          AlignedVector<Number>                 &tmp,
          const unsigned int                     n_rows_1d_non_templated,
          const std::array<const Number *, dim> &mass_matrix,
          const std::array<const Number *, dim> &derivative_matrix)
    {
      const unsigned int n_rows_1d = n_rows_1d_templated == 0 ?
                                       n_rows_1d_non_templated :
                                       n_rows_1d_templated;
      const unsigned int n         = Utilities::fixed_power<dim>(n_rows_1d);

      tmp.resize_fast(n * 2);
      Number *t = tmp.begin();

      internal::EvaluatorTensorProduct<internal::evaluate_general,
                                       dim,
                                       n_rows_1d_templated,
                                       n_rows_1d_templated,
                                       Number>
        eval({}, {}, {}, n_rows_1d, n_rows_1d);

      if (dim == 1)
        {
          const Number *A = derivative_matrix[0];
          eval.template apply<0, false, false>(A, src, dst);
        }

      else if (dim == 2)
        {
          const Number *A0 = derivative_matrix[0];
          const Number *M0 = mass_matrix[0];
          const Number *A1 = derivative_matrix[1];
          const Number *M1 = mass_matrix[1];
          eval.template apply<0, false, false>(M0, src, t);
          eval.template apply<1, false, false>(A1, t, dst);
          eval.template apply<0, false, false>(A0, src, t);
          eval.template apply<1, false, true>(M1, t, dst);
        }

      else if (dim == 3)
        {
          const Number *A0 = derivative_matrix[0];
          const Number *M0 = mass_matrix[0];
          const Number *A1 = derivative_matrix[1];
          const Number *M1 = mass_matrix[1];
          const Number *A2 = derivative_matrix[2];
          const Number *M2 = mass_matrix[2];
          eval.template apply<0, false, false>(M0, src, t + n);
          eval.template apply<1, false, false>(M1, t + n, t);
          eval.template apply<2, false, false>(A2, t, dst);
          eval.template apply<1, false, false>(A1, t + n, t);
          eval.template apply<0, false, false>(A0, src, t + n);
          eval.template apply<1, false, true>(M1, t + n, t);
          eval.template apply<2, false, true>(M2, t, dst);
        }

      else
        AssertThrow(false, ExcNotImplemented());
    }



    template <int n_rows_1d_templated, std::size_t dim, typename Number>
    void
    apply_inverse(Number            *dst,
                  const Number      *src,
                  const unsigned int n_rows_1d_non_templated,
                  const std::array<const Number *, dim> &eigenvectors,
                  const std::array<const Number *, dim> &eigenvalues,
                  const Number *inverted_eigenvalues = nullptr)
    {
      const unsigned int n_rows_1d = n_rows_1d_templated == 0 ?
                                       n_rows_1d_non_templated :
                                       n_rows_1d_templated;

      internal::EvaluatorTensorProduct<internal::evaluate_general,
                                       dim,
                                       n_rows_1d_templated,
                                       n_rows_1d_templated,
                                       Number>
        eval({}, {}, {}, n_rows_1d, n_rows_1d);

      // NOTE: dof_to_quad has to be interpreted as 'dof to eigenvalue index'
      //       --> apply<.,true,.> (S,src,dst) calculates dst = S^T * src,
      //       --> apply<.,false,.> (S,src,dst) calculates dst = S * src,
      //       while the eigenvectors are stored column-wise in S, i.e.
      //       rows correspond to dofs whereas columns to eigenvalue indices!
      if (dim == 1)
        {
          const Number *S = eigenvectors[0];
          eval.template apply<0, true, false>(S, src, dst);

          for (unsigned int i = 0; i < n_rows_1d; ++i)
            if (inverted_eigenvalues)
              dst[i] *= inverted_eigenvalues[i];
            else
              dst[i] /= eigenvalues[0][i];

          eval.template apply<0, false, false>(S, dst, dst);
        }

      else if (dim == 2)
        {
          const Number *S0 = eigenvectors[0];
          const Number *S1 = eigenvectors[1];
          eval.template apply<0, true, false>(S0, src, dst);
          eval.template apply<1, true, false>(S1, dst, dst);

          for (unsigned int i1 = 0, c = 0; i1 < n_rows_1d; ++i1)
            for (unsigned int i0 = 0; i0 < n_rows_1d; ++i0, ++c)
              if (inverted_eigenvalues)
                dst[c] *= inverted_eigenvalues[c];
              else
                dst[c] /= (eigenvalues[1][i1] + eigenvalues[0][i0]);

          eval.template apply<1, false, false>(S1, dst, dst);
          eval.template apply<0, false, false>(S0, dst, dst);
        }

      else if (dim == 3)
        {
          const Number *S0 = eigenvectors[0];
          const Number *S1 = eigenvectors[1];
          const Number *S2 = eigenvectors[2];
          eval.template apply<0, true, false>(S0, src, dst);
          eval.template apply<1, true, false>(S1, dst, dst);
          eval.template apply<2, true, false>(S2, dst, dst);

          for (unsigned int i2 = 0, c = 0; i2 < n_rows_1d; ++i2)
            for (unsigned int i1 = 0; i1 < n_rows_1d; ++i1)
              for (unsigned int i0 = 0; i0 < n_rows_1d; ++i0, ++c)
                if (inverted_eigenvalues)
                  dst[c] *= inverted_eigenvalues[c];
                else
                  dst[c] /= (eigenvalues[2][i2] + eigenvalues[1][i1] +
                             eigenvalues[0][i0]);

          eval.template apply<2, false, false>(S2, dst, dst);
          eval.template apply<1, false, false>(S1, dst, dst);
          eval.template apply<0, false, false>(S0, dst, dst);
        }

      else
        DEAL_II_NOT_IMPLEMENTED();
    }



    template <int n_rows_1d_templated, std::size_t dim, typename Number>
    void
    select_vmult(Number                                *dst,
                 const Number                          *src,
                 AlignedVector<Number>                 &tmp,
                 const unsigned int                     n_rows_1d,
                 const std::array<const Number *, dim> &mass_matrix,
                 const std::array<const Number *, dim> &derivative_matrix);



    template <int n_rows_1d_templated, std::size_t dim, typename Number>
    void
    select_apply_inverse(Number                                *dst,
                         const Number                          *src,
                         const unsigned int                     n_rows_1d,
                         const std::array<const Number *, dim> &eigenvectors,
                         const std::array<const Number *, dim> &eigenvalues,
                         const Number *inverted_eigenvalues = nullptr);
  } // namespace TensorProductMatrixSymmetricSum
} // namespace internal


template <int dim, typename Number, int n_rows_1d>
inline unsigned int
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::m() const
{
  unsigned int m = mass_matrix[0].n_rows();
  for (unsigned int d = 1; d < dim; ++d)
    m *= mass_matrix[d].n_rows();
  return m;
}



template <int dim, typename Number, int n_rows_1d>
inline unsigned int
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::n() const
{
  unsigned int n = mass_matrix[0].n_cols();
  for (unsigned int d = 1; d < dim; ++d)
    n *= mass_matrix[d].n_cols();
  return n;
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::vmult(
  const ArrayView<Number>       &dst_view,
  const ArrayView<const Number> &src_view) const
{
  std::lock_guard<std::mutex> lock(this->mutex);
  this->vmult(dst_view, src_view, this->tmp_array);
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::vmult(
  const ArrayView<Number>       &dst_view,
  const ArrayView<const Number> &src_view,
  AlignedVector<Number>         &tmp_array) const
{
  AssertDimension(dst_view.size(), this->m());
  AssertDimension(src_view.size(), this->n());

  Number       *dst = dst_view.begin();
  const Number *src = src_view.begin();

  std::array<const Number *, dim> mass_matrix, derivative_matrix;

  for (unsigned int d = 0; d < dim; ++d)
    {
      mass_matrix[d]       = &this->mass_matrix[d](0, 0);
      derivative_matrix[d] = &this->derivative_matrix[d](0, 0);
    }

  const unsigned int n_rows_1d_non_templated = this->mass_matrix[0].n_rows();

  if (n_rows_1d != -1)
    internal::TensorProductMatrixSymmetricSum::vmult<
      n_rows_1d == -1 ? 0 : n_rows_1d>(dst,
                                       src,
                                       tmp_array,
                                       n_rows_1d_non_templated,
                                       mass_matrix,
                                       derivative_matrix);
  else
    internal::TensorProductMatrixSymmetricSum::select_vmult<1>(
      dst,
      src,
      tmp_array,
      n_rows_1d_non_templated,
      mass_matrix,
      derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::apply_inverse(
  const ArrayView<Number>       &dst_view,
  const ArrayView<const Number> &src_view) const
{
  AssertDimension(dst_view.size(), this->n());
  AssertDimension(src_view.size(), this->m());

  Number       *dst = dst_view.begin();
  const Number *src = src_view.begin();

  std::array<const Number *, dim> eigenvectors, eigenvalues;

  for (unsigned int d = 0; d < dim; ++d)
    {
      eigenvectors[d] = &this->eigenvectors[d](0, 0);
      eigenvalues[d]  = this->eigenvalues[d].data();
    }

  const unsigned int n_rows_1d_non_templated = this->mass_matrix[0].n_rows();

  if (n_rows_1d != -1)
    internal::TensorProductMatrixSymmetricSum::apply_inverse<
      n_rows_1d == -1 ? 0 : n_rows_1d>(
      dst, src, n_rows_1d_non_templated, eigenvectors, eigenvalues);
  else
    internal::TensorProductMatrixSymmetricSum::select_apply_inverse<1>(
      dst, src, n_rows_1d_non_templated, eigenvectors, eigenvalues);
}



template <int dim, typename Number, int n_rows_1d>
std::size_t
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::memory_consumption()
  const
{
  return MemoryConsumption::memory_consumption(mass_matrix) +
         MemoryConsumption::memory_consumption(derivative_matrix) +
         MemoryConsumption::memory_consumption(eigenvalues) +
         MemoryConsumption::memory_consumption(eigenvectors) +
         MemoryConsumption::memory_consumption(tmp_array);
}



template <int dim, typename Number, int n_rows_1d>
template <typename T>
inline TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::
  TensorProductMatrixSymmetricSum(const T &mass_matrix,
                                  const T &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
template <typename T>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::reinit(
  const T &mass_matrix,
  const T &derivative_matrix)
{
  this->mass_matrix =
    internal::TensorProductMatrixSymmetricSum::convert<dim>(mass_matrix);
  this->derivative_matrix =
    internal::TensorProductMatrixSymmetricSum::convert<dim>(derivative_matrix);

  internal::TensorProductMatrixSymmetricSum::setup(this->mass_matrix,
                                                   this->derivative_matrix,
                                                   this->eigenvectors,
                                                   this->eigenvalues);
}



template <int dim, typename Number, int n_rows_1d>
TensorProductMatrixSymmetricSumCollection<dim, Number, n_rows_1d>::
  AdditionalData::AdditionalData(const bool compress_matrices,
                                 const bool precompute_inverse_diagonal)
  : compress_matrices(compress_matrices)
  , precompute_inverse_diagonal(precompute_inverse_diagonal)
{}



template <int dim, typename Number, int n_rows_1d>
TensorProductMatrixSymmetricSumCollection<dim, Number, n_rows_1d>::
  TensorProductMatrixSymmetricSumCollection(
    const AdditionalData &additional_data)
  : compress_matrices(additional_data.compress_matrices)
  , precompute_inverse_diagonal(additional_data.precompute_inverse_diagonal)
{}



template <int dim, typename Number, int n_rows_1d>
void
TensorProductMatrixSymmetricSumCollection<dim, Number, n_rows_1d>::reserve(
  const unsigned int size)
{
  if (compress_matrices == false)
    mass_and_derivative_matrices.resize(size * dim);
  else
    indices.assign(size * dim, numbers::invalid_unsigned_int);
}



template <int dim, typename Number, int n_rows_1d>
template <typename T>
void
TensorProductMatrixSymmetricSumCollection<dim, Number, n_rows_1d>::insert(
  const unsigned int index,
  const T           &Ms_in,
  const T           &Ks_in)
{
  const auto Ms =
    internal::TensorProductMatrixSymmetricSum::convert<dim>(Ms_in);
  const auto Ks =
    internal::TensorProductMatrixSymmetricSum::convert<dim>(Ks_in);

  for (unsigned int d = 0; d < dim; ++d)
    {
      if (compress_matrices == false)
        {
          const MatrixPairType matrix(Ms[d], Ks[d]);
          mass_and_derivative_matrices[index * dim + d] = matrix;
        }
      else
        {
          using VectorizedArrayTrait =
            dealii::internal::VectorizedArrayTrait<Number>;

          std::bitset<VectorizedArrayTrait::width()> mask;

          for (unsigned int v = 0; v < VectorizedArrayTrait::width(); ++v)
            {
              typename VectorizedArrayTrait::value_type a = 0.0;

              for (unsigned int i = 0; i < Ms[d].size(0); ++i)
                for (unsigned int j = 0; j < Ms[d].size(1); ++j)
                  {
                    a += std::abs(VectorizedArrayTrait::get(Ms[d][i][j], v));
                    a += std::abs(VectorizedArrayTrait::get(Ks[d][i][j], v));
                  }

              mask[v] = (a != 0.0);
            }

          const MatrixPairTypeWithMask matrix{mask, {Ms[d], Ks[d]}};

          const auto ptr = cache.find(matrix);

          if (ptr != cache.end())
            {
              const auto ptr_index     = ptr->second;
              indices[index * dim + d] = ptr_index;

              if ([&]() {
                    for (unsigned int v = 0; v < VectorizedArrayTrait::width();
                         ++v)
                      if ((mask[v] == true) && (ptr->first.first[v] == false))
                        return false;

                    return true;
                  }())
                {
                  // nothing to do
                }
              else
                {
                  auto mask_new = ptr->first.first;
                  auto Ms_new   = ptr->first.second.first;
                  auto Ks_new   = ptr->first.second.second;

                  for (unsigned int v = 0; v < VectorizedArrayTrait::width();
                       ++v)
                    if (mask_new[v] == false && mask[v] == true)
                      {
                        mask_new[v] = true;

                        for (unsigned int i = 0; i < Ms_new.size(0); ++i)
                          for (unsigned int j = 0; j < Ms_new.size(1); ++j)
                            {
                              VectorizedArrayTrait::get(Ms_new[i][j], v) =
                                VectorizedArrayTrait::get(Ms[d][i][j], v);
                              VectorizedArrayTrait::get(Ks_new[i][j], v) =
                                VectorizedArrayTrait::get(Ks[d][i][j], v);
                            }
                      }

                  cache.erase(ptr);

                  const MatrixPairTypeWithMask entry_new{mask_new,
                                                         {Ms_new, Ks_new}};

                  const auto ptr_ = cache.find(entry_new);
                  AssertThrow(ptr_ == cache.end(), ExcNotImplemented());

                  cache[entry_new] = ptr_index;
                }
            }
          else
            {
              const auto size          = cache.size();
              indices[index * dim + d] = size;
              cache[matrix]            = size;
            }
        }
    }
}



template <int dim, typename Number, int n_rows_1d>
void
TensorProductMatrixSymmetricSumCollection<dim, Number, n_rows_1d>::finalize()
{
  const auto store = [&](const unsigned int    index,
                         const MatrixPairType &M_and_K) {
    std::array<Table<2, Number>, 1> mass_matrix;
    mass_matrix[0] = M_and_K.first;

    std::array<Table<2, Number>, 1> derivative_matrix;
    derivative_matrix[0] = M_and_K.second;

    std::array<Table<2, Number>, 1>      eigenvectors;
    std::array<AlignedVector<Number>, 1> eigenvalues;

    internal::TensorProductMatrixSymmetricSum::setup(mass_matrix,
                                                     derivative_matrix,
                                                     eigenvectors,
                                                     eigenvalues);

    for (unsigned int i = 0, m = matrix_ptr[index], v = vector_ptr[index];
         i < mass_matrix[0].n_rows();
         ++i, ++v)
      {
        for (unsigned int j = 0; j < mass_matrix[0].n_cols(); ++j, ++m)
          {
            this->mass_matrices[m]       = mass_matrix[0][i][j];
            this->derivative_matrices[m] = derivative_matrix[0][i][j];
            this->eigenvectors[m]        = eigenvectors[0][i][j];
          }

        this->eigenvalues[v] = eigenvalues[0][i];
      }
  };

  if (compress_matrices == false)
    {
      // case 1) no compression requested

      AssertDimension(cache.size(), 0);
      AssertDimension(indices.size(), 0);

      this->vector_ptr.resize(mass_and_derivative_matrices.size() + 1);
      this->matrix_ptr.resize(mass_and_derivative_matrices.size() + 1);

      for (unsigned int i = 0; i < mass_and_derivative_matrices.size(); ++i)
        {
          const auto &M = mass_and_derivative_matrices[i].first;

          this->vector_ptr[i + 1] = M.n_rows();
          this->matrix_ptr[i + 1] = M.n_rows() * M.n_cols();
        }

      for (unsigned int i = 0; i < mass_and_derivative_matrices.size(); ++i)
        {
          this->vector_ptr[i + 1] += this->vector_ptr[i];
          this->matrix_ptr[i + 1] += this->matrix_ptr[i];
        }

      this->mass_matrices.resize_fast(matrix_ptr.back());
      this->derivative_matrices.resize_fast(matrix_ptr.back());
      this->eigenvectors.resize_fast(matrix_ptr.back());
      this->eigenvalues.resize_fast(vector_ptr.back());

      for (unsigned int i = 0; i < mass_and_derivative_matrices.size(); ++i)
        store(i, mass_and_derivative_matrices[i]);

      mass_and_derivative_matrices.clear();
    }
  else if (cache.size() == indices.size())
    {
      // case 2) compression requested but none possible

      this->vector_ptr.resize(cache.size() + 1);
      this->matrix_ptr.resize(cache.size() + 1);

      std::map<unsigned int, MatrixPairType> inverted_cache;

      for (const auto &i : cache)
        inverted_cache[i.second] = i.first.second;

      for (unsigned int i = 0; i < indices.size(); ++i)
        {
          const auto &M = inverted_cache[indices[i]].first;

          this->vector_ptr[i + 1] = M.n_rows();
          this->matrix_ptr[i + 1] = M.n_rows() * M.n_cols();
        }

      for (unsigned int i = 0; i < cache.size(); ++i)
        {
          this->vector_ptr[i + 1] += this->vector_ptr[i];
          this->matrix_ptr[i + 1] += this->matrix_ptr[i];
        }

      this->mass_matrices.resize_fast(matrix_ptr.back());
      this->derivative_matrices.resize_fast(matrix_ptr.back());
      this->eigenvectors.resize_fast(matrix_ptr.back());
      this->eigenvalues.resize_fast(vector_ptr.back());

      for (unsigned int i = 0; i < indices.size(); ++i)
        store(i, inverted_cache[indices[i]]);

      indices.clear();
      cache.clear();
    }
  else
    {
      // case 3) compress

      this->vector_ptr.resize(cache.size() + 1);
      this->matrix_ptr.resize(cache.size() + 1);

      for (const auto &i : cache)
        {
          const auto &M = i.first.second.first;

          this->vector_ptr[i.second + 1] = M.n_rows();
          this->matrix_ptr[i.second + 1] = M.n_rows() * M.n_cols();
        }

      for (unsigned int i = 0; i < cache.size(); ++i)
        {
          this->vector_ptr[i + 1] += this->vector_ptr[i];
          this->matrix_ptr[i + 1] += this->matrix_ptr[i];
        }

      this->mass_matrices.resize_fast(matrix_ptr.back());
      this->derivative_matrices.resize_fast(matrix_ptr.back());
      this->eigenvectors.resize_fast(matrix_ptr.back());
      this->eigenvalues.resize_fast(vector_ptr.back());

      for (const auto &i : cache)
        store(i.second, i.first.second);

      cache.clear();
    }

  if (precompute_inverse_diagonal)
    {
      if (dim == 1)
        {
          // 1D case: simply invert 1D eigenvalues
          for (unsigned int i = 0; i < this->eigenvalues.size(); ++i)
            this->eigenvalues[i] = Number(1.0) / this->eigenvalues[i];
          std::swap(this->inverted_eigenvalues, eigenvalues);
        }
      else
        {
          // 2D and 3D case: we have 2 or 3 1d eigenvalues so that we
          // need to combine these

          // step 1) if eigenvalues/eigenvectors are compressed, we
          // need to compress the diagonal (the combination of ev
          // indices) as well. This is an optional step.
          std::vector<unsigned int> indices_ev;

          if (indices.size() > 0)
            {
              // 1a) create cache (ev indics -> diag index)
              const unsigned int n_cells = indices.size() / dim;
              std::map<std::array<unsigned int, dim>, unsigned int> cache_ev;
              std::vector<unsigned int> cache_ev_idx(n_cells);

              for (unsigned int i = 0, c = 0; i < n_cells; ++i)
                {
                  std::array<unsigned int, dim> id;

                  for (unsigned int d = 0; d < dim; ++d, ++c)
                    id[d] = indices[c];

                  const auto id_ptr = cache_ev.find(id);

                  if (id_ptr == cache_ev.end())
                    {
                      const auto size = cache_ev.size();
                      cache_ev_idx[i] = size;
                      cache_ev[id]    = size;
                    }
                  else
                    {
                      cache_ev_idx[i] = id_ptr->second;
                    }
                }

              // 1b) store diagonal indices for each cell
              std::vector<unsigned int> new_indices;
              new_indices.reserve(indices.size() / dim * (dim + 1));

              for (unsigned int i = 0, c = 0; i < n_cells; ++i)
                {
                  for (unsigned int d = 0; d < dim; ++d, ++c)
                    new_indices.push_back(indices[c]);
                  new_indices.push_back(cache_ev_idx[i]);
                }

              // 1c) transpose cache (diag index -> ev indices)
              indices_ev.resize(cache_ev.size() * dim);
              for (const auto &entry : cache_ev)
                for (unsigned int d = 0; d < dim; ++d)
                  indices_ev[entry.second * dim + d] = entry.first[d];

              std::swap(this->indices, new_indices);
            }

          // step 2) allocate memory and set pointers
          const unsigned int n_diag =
            ((indices_ev.size() > 0) ? indices_ev.size() :
                                       (matrix_ptr.size() - 1)) /
            dim;

          std::vector<unsigned int> new_vector_ptr(n_diag + 1, 0);
          std::vector<unsigned int> new_vector_n_rows_1d(n_diag, 0);

          for (unsigned int i = 0; i < n_diag; ++i)
            {
              const unsigned int c = (indices_ev.size() > 0) ?
                                       indices_ev[dim * i + 0] :
                                       (dim * i + 0);

              const unsigned int n_rows = vector_ptr[c + 1] - vector_ptr[c];

              new_vector_n_rows_1d[i] = n_rows;
              new_vector_ptr[i + 1]   = Utilities::pow(n_rows, dim);
            }

          for (unsigned int i = 0; i < n_diag; ++i)
            new_vector_ptr[i + 1] += new_vector_ptr[i];

          this->inverted_eigenvalues.resize(new_vector_ptr.back());

          // step 3) loop over all unique diagonal entries and invert
          for (unsigned int i = 0; i < n_diag; ++i)
            {
              std::array<Number *, dim> evs;

              for (unsigned int d = 0; d < dim; ++d)
                evs[d] =
                  &this
                     ->eigenvalues[this->vector_ptr[(indices_ev.size() > 0) ?
                                                      indices_ev[dim * i + d] :
                                                      (dim * i + d)]];

              const unsigned int mm = new_vector_n_rows_1d[i];
              if (dim == 2)
                {
                  for (unsigned int i1 = 0, c = 0; i1 < mm; ++i1)
                    for (unsigned int i0 = 0; i0 < mm; ++i0, ++c)
                      this->inverted_eigenvalues[new_vector_ptr[i] + c] =
                        Number(1.0) / (evs[1][i1] + evs[0][i0]);
                }
              else
                {
                  for (unsigned int i2 = 0, c = 0; i2 < mm; ++i2)
                    for (unsigned int i1 = 0; i1 < mm; ++i1)
                      for (unsigned int i0 = 0; i0 < mm; ++i0, ++c)
                        this->inverted_eigenvalues[new_vector_ptr[i] + c] =
                          Number(1.0) / (evs[2][i2] + evs[1][i1] + evs[0][i0]);
                }
            }

          // step 4) clean up
          std::swap(this->vector_ptr, new_vector_ptr);
          std::swap(this->vector_n_rows_1d, new_vector_n_rows_1d);
        }

      this->eigenvalues.clear();
    }
}



template <int dim, typename Number, int n_rows_1d>
void
TensorProductMatrixSymmetricSumCollection<dim, Number, n_rows_1d>::
  apply_inverse(const unsigned int             index,
                const ArrayView<Number>       &dst_in,
                const ArrayView<const Number> &src_in) const
{
  Number       *dst = dst_in.begin();
  const Number *src = src_in.begin();

  if (this->eigenvalues.empty() == false)
    {
      std::array<const Number *, dim> eigenvectors;
      std::array<const Number *, dim> eigenvalues;
      unsigned int                    n_rows_1d_non_templated = 0;

      for (unsigned int d = 0; d < dim; ++d)
        {
          const unsigned int translated_index =
            (indices.size() > 0) ? indices[dim * index + d] : (dim * index + d);

          eigenvectors[d] =
            this->eigenvectors.data() + matrix_ptr[translated_index];
          eigenvalues[d] =
            this->eigenvalues.data() + vector_ptr[translated_index];
          n_rows_1d_non_templated =
            vector_ptr[translated_index + 1] - vector_ptr[translated_index];
        }

      if (n_rows_1d != -1)
        internal::TensorProductMatrixSymmetricSum::apply_inverse<
          n_rows_1d == -1 ? 0 : n_rows_1d>(
          dst, src, n_rows_1d_non_templated, eigenvectors, eigenvalues);
      else
        internal::TensorProductMatrixSymmetricSum::select_apply_inverse<1>(
          dst, src, n_rows_1d_non_templated, eigenvectors, eigenvalues);
    }
  else
    {
      std::array<const Number *, dim> eigenvectors;
      const Number                   *inverted_eigenvalues    = nullptr;
      unsigned int                    n_rows_1d_non_templated = 0;

      for (unsigned int d = 0; d < dim; ++d)
        {
          const unsigned int translated_index =
            (indices.size() > 0) ?
              indices[((dim == 1) ? 1 : (dim + 1)) * index + d] :
              (dim * index + d);

          eigenvectors[d] =
            this->eigenvectors.data() + matrix_ptr[translated_index];
        }

      {
        const unsigned int translated_index =
          ((indices.size() > 0) && (dim != 1)) ?
            indices[(dim + 1) * index + dim] :
            index;

        inverted_eigenvalues =
          this->inverted_eigenvalues.data() + vector_ptr[translated_index];
        n_rows_1d_non_templated =
          (dim == 1) ?
            (vector_ptr[translated_index + 1] - vector_ptr[translated_index]) :
            vector_n_rows_1d[translated_index];
      }

      if (n_rows_1d != -1)
        internal::TensorProductMatrixSymmetricSum::apply_inverse<
          n_rows_1d == -1 ? 0 : n_rows_1d>(dst,
                                           src,
                                           n_rows_1d_non_templated,
                                           eigenvectors,
                                           {},
                                           inverted_eigenvalues);
      else
        internal::TensorProductMatrixSymmetricSum::select_apply_inverse<1>(
          dst,
          src,
          n_rows_1d_non_templated,
          eigenvectors,
          {},
          inverted_eigenvalues);
    }
}



template <int dim, typename Number, int n_rows_1d>
std::size_t
TensorProductMatrixSymmetricSumCollection<dim, Number, n_rows_1d>::
  memory_consumption() const
{
  return MemoryConsumption::memory_consumption(indices) +
         MemoryConsumption::memory_consumption(mass_matrices) +
         MemoryConsumption::memory_consumption(derivative_matrices) +
         MemoryConsumption::memory_consumption(eigenvectors) +
         MemoryConsumption::memory_consumption(eigenvalues) +
         MemoryConsumption::memory_consumption(matrix_ptr) +
         MemoryConsumption::memory_consumption(vector_ptr);
}



template <int dim, typename Number, int n_rows_1d>
std::size_t
TensorProductMatrixSymmetricSumCollection<dim, Number, n_rows_1d>::
  storage_size() const
{
  if (matrix_ptr.empty())
    return 0; // if not initialized

  return matrix_ptr.size() - 1;
}



#endif

DEAL_II_NAMESPACE_CLOSE

#endif
