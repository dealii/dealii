// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_physics_notation_h
#define dealii_physics_notation_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <type_traits>

DEAL_II_NAMESPACE_OPEN


namespace Physics
{
  namespace Notation
  {
    /**
     * @brief A namespace with functions that assist in the conversion of
     * vectors and tensors to and from a compressed format using Kelvin notation
     * and weighting.
     *
     * Both Kelvin and Voigt notation adopt the same indexing convention.
     * With specific reference to the spatial dimension 3 case, for
     * a rank-2 symmetric tensor $\mathbf{S}$ we enumerate its tensor
     * components
     * @f[
     * \mathbf{S} \dealcoloneq \left[ \begin{array}{ccc}
     *  S_{00}          & S_{01}          & S_{02} \\
     *  S_{10} = S_{01} & S_{11}          & S_{12} \\
     *  S_{20} = S_{02} & S_{21} = S_{12} & S_{22}
     * \end{array} \right]
     * \quad \Rightarrow \quad
     * \left[ \begin{array}{ccc}
     *  n = 0 & n = 5 & n = 4 \\
     *  sym   & n = 1 & n = 3 \\
     *  sym   & sym   & n = 2
     * \end{array} \right] ,
     * @f]
     * where $n$ denotes the Kelvin index for the tensor component,
     * while for a general rank-2 tensor $\mathbf{T}$
     * @f[
     * \mathbf{T} \dealcoloneq \left[ \begin{array}{ccc}
     *  T_{00} & T_{01} & T_{02} \\
     *  T_{10} & T_{11} & T_{12} \\
     *  T_{20} & T_{21} & T_{22}
     * \end{array}\right]
     * \quad \Rightarrow \quad
     * \left[ \begin{array}{ccc}
     *  n = 0 & n = 5 & n = 4 \\
     *  n = 6 & n = 1 & n = 3 \\
     *  n = 7 & n = 8 & n = 2
     * \end{array}\right] ,
     * @f]
     * and for a rank-1 tensor $\mathbf{v}$
     * @f[
     * \mathbf{v} \dealcoloneq \left[ \begin{array}{c}
     *  v_{0} \\ v_{1} \\ v_{2}
     * \end{array}\right]
     * \quad \Rightarrow \quad
     * \left[ \begin{array}{c}
     *  n = 0 \\ n = 1 \\ n = 2
     * \end{array}\right] .
     * @f]
     * To summarize, the relationship between tensor and Kelvin indices for both
     * the three-dimensional case and the analogously discerned two-dimensional
     * case outlined in the following table:
     * <table>
     * <tr>
     *   <th align="center"> Dimension 2 </th>
     *   <th align="center"> Dimension 3 </th>
     * </tr>
     * <tr>
     * <td align="middle">
     *   <table>
     *   <tr>
     *     <th>Tensor index pairs</th>
     *     <th>Kelvin index</th>
     *   </tr>
     *   <tr>
     *     <td align="center">00</td>
     *     <td align="center">0</td>
     *   </tr>
     *   <tr>
     *     <td align="center">11</td>
     *     <td align="center">1</td>
     *   </tr>
     *   <tr>
     *     <td align="center">12</td>
     *     <td align="center">2</td>
     *   </tr>
     *   <tr>
     *     <td align="center">21</td>
     *     <td align="center">3</td>
     *   </tr>
     *   </table>
     * </td>
     * <td align="middle">
     *   <table>
     *   <tr>
     *     <th>Tensor index pairs</th>
     *     <th>Kelvin index</th>
     *   </tr>
     *   <tr>
     *     <td align="center">00</td>
     *     <td align="center">0</td>
     *   </tr>
     *   <tr>
     *     <td align="center">11</td>
     *     <td align="center">1</td>
     *   </tr>
     *   <tr>
     *     <td align="center">22</td>
     *     <td align="center">2</td>
     *   </tr>
     *   <tr>
     *     <td align="center">12</td>
     *     <td align="center">3</td>
     *   </tr>
     *   <tr>
     *     <td align="center">02</td>
     *     <td align="center">4</td>
     *   </tr>
     *   <tr>
     *     <td align="center">01</td>
     *     <td align="center">5</td>
     *   </tr>
     *   <tr>
     *     <td align="center">10</td>
     *     <td align="center">6</td>
     *   </tr>
     *   <tr>
     *     <td align="center">20</td>
     *     <td align="center">7</td>
     *   </tr>
     *   <tr>
     *     <td align="center">21</td>
     *     <td align="center">8</td>
     *   </tr>
     *   </table>
     * </td>
     * </tr>
     * </table>
     *
     * To illustrate the purpose of this notation, consider the rank-2 symmetric
     * tensors $\mathbf{S}$ and $\mathbf{E}$ that are related to one another by
     * $\mathbf{S} = \cal{C} : \mathbf{E}$, where the operator $\cal{C}$ is a
     * fourth-order symmetric tensor. As opposed to the commonly used Voigt
     * notation, Kelvin (or Mandel) notation keeps the same definition of the
     * inner product $\mathbf{S} : \mathbf{E}$ when both $\mathbf{S}$ and
     * $\mathbf{E}$ are symmetric. In general, the inner product of all
     * symmetric and general tensors remain the same regardless of the notation
     * with which it is represented.
     *
     * To achieve these two properties, namely that
     * @f[
     * \mathbf{S} = \cal{C} : \mathbf{E}
     * \quad \Rightarrow   \quad
     * \tilde{\mathbf{S}} = \tilde{\cal{C}} \; \tilde{\mathbf{E}}
     * @f]
     * and
     * @f[
     * \mathbf{S} : \mathbf{E}
     * \, \equiv \,
     * \tilde{\mathbf{S}} \cdot \tilde{\mathbf{E}} ,
     * @f]
     * it holds that the Kelvin-condensed equivalents of the previously defined
     * symmetric tensors, indicated by the $\tilde{\left(\bullet\right)}$, must
     * be defined as
     * @f[
     * \tilde{\mathbf{S}}
     *   = \left[ \begin{array}{c}
     *   S_{00} \\ S_{11} \\ S_{22} \\ \sqrt{2} S_{12} \\ \sqrt{2} S_{02} \\
     * \sqrt{2} S_{01} \end{array}\right] \quad \text{and} \quad
     * \tilde{\mathbf{E}}
     *   = \left[ \begin{array}{c}
     *   E_{00} \\ E_{11} \\ E_{22} \\ \sqrt{2} E_{12} \\ \sqrt{2} E_{02} \\
     * \sqrt{2} E_{01} \end{array}\right] .
     * @f]
     * The corresponding and consistent condensed fourth-order symmetric tensor
     * is
     * @f[
     * \tilde{\cal{C}}
     *   = \left[ \begin{array}{cccccc}
     *   \tilde{\cal{C}}_{00} & \tilde{\cal{C}}_{01} & \tilde{\cal{C}}_{02} &
     * \tilde{\cal{C}}_{03} & \tilde{\cal{C}}_{04} & \tilde{\cal{C}}_{05} \\
     *   \tilde{\cal{C}}_{10} & \tilde{\cal{C}}_{11} & \tilde{\cal{C}}_{12} &
     * \tilde{\cal{C}}_{13} & \tilde{\cal{C}}_{14} & \tilde{\cal{C}}_{15} \\
     *   \tilde{\cal{C}}_{20} & \tilde{\cal{C}}_{21} & \tilde{\cal{C}}_{22} &
     * \tilde{\cal{C}}_{23} & \tilde{\cal{C}}_{24} & \tilde{\cal{C}}_{25} \\
     *   \tilde{\cal{C}}_{30} & \tilde{\cal{C}}_{31} & \tilde{\cal{C}}_{32} &
     * \tilde{\cal{C}}_{33} & \tilde{\cal{C}}_{34} & \tilde{\cal{C}}_{35} \\
     *   \tilde{\cal{C}}_{40} & \tilde{\cal{C}}_{41} & \tilde{\cal{C}}_{42} &
     * \tilde{\cal{C}}_{43} & \tilde{\cal{C}}_{44} & \tilde{\cal{C}}_{45} \\
     *   \tilde{\cal{C}}_{50} & \tilde{\cal{C}}_{51} & \tilde{\cal{C}}_{52} &
     * \tilde{\cal{C}}_{53} & \tilde{\cal{C}}_{54} & \tilde{\cal{C}}_{55}
     *   \end{array}\right]
     *   \equiv
     *   \left[ \begin{array}{cccccc}
     *   {\cal{C}}_{0000}           & {\cal{C}}_{0011}          &
     * {\cal{C}}_{0022}           & \sqrt{2} {\cal{C}}_{0012}  & \sqrt{2}
     * {\cal{C}}_{0002}  & \sqrt{2} {\cal{C}}_{0001} \\
     *   {\cal{C}}_{1100}           & {\cal{C}}_{1111}          &
     * {\cal{C}}_{1122}           & \sqrt{2} {\cal{C}}_{1112}  & \sqrt{2}
     * {\cal{C}}_{1102}  & \sqrt{2} {\cal{C}}_{1101} \\
     *   {\cal{C}}_{2200}           & {\cal{C}}_{2211}          &
     * {\cal{C}}_{2222}           & \sqrt{2} {\cal{C}}_{2212}  & \sqrt{2}
     * {\cal{C}}_{2202}  & \sqrt{2} {\cal{C}}_{2201} \\
     *   \sqrt{2} {\cal{C}}_{1200}  & \sqrt{2} {\cal{C}}_{1211} & \sqrt{2}
     * {\cal{C}}_{1222}  & 2 {\cal{C}}_{1212}         & 2 {\cal{C}}_{1202} & 2
     * {\cal{C}}_{1201}        \\
     *   \sqrt{2} {\cal{C}}_{0200}  & \sqrt{2} {\cal{C}}_{0211} & \sqrt{2}
     * {\cal{C}}_{0222}  & 2 {\cal{C}}_{0212}         & 2 {\cal{C}}_{0202} & 2
     * {\cal{C}}_{0201}        \\ \sqrt{2} {\cal{C}}_{0100}  & \sqrt{2}
     * {\cal{C}}_{0111} & \sqrt{2} {\cal{C}}_{0122}  & 2 {\cal{C}}_{0112} & 2
     * {\cal{C}}_{0102}         & 2 {\cal{C}}_{0101} \end{array}\right] .
     * @f]
     * The mapping from the two Kelvin indices of the FullMatrix
     * $\tilde{\cal{C}}$ to the rank-4 SymmetricTensor $\cal{C}$ can be inferred
     * using the table shown above.
     *
     * An important observation is that both the left-hand side tensor
     * $\tilde{\mathbf{S}}$ and right-hand side tensor $\tilde{\mathbf{E}}$ have
     * the same form; this is a property that is not present in Voigt notation.
     * The various factors introduced into $\tilde{\mathbf{S}}$,
     * $\tilde{\mathbf{E}}$ and $\tilde{\cal{C}}$ account for the symmetry of
     * the tensors. The Kelvin description of their non-symmetric counterparts
     * include no such factors.
     *
     * Some useful references that show how this notation works include, amongst
     * others,
     * @code{.bib}
     * @article{Nagel2016,
     *   author  = {Nagel, T. and G{\"o}rke, U-J. and Moerman, K. and Kolditz,
     *              O.},
     *   title   = {On advantages of the Kelvin mapping in finite element
     *              implementations of deformation processes},
     *   journal = {Environmental Earth Sciences},
     *   year    = {2016},
     *   volume  = {75},
     *   number  = {11},
     *   pages   = {937}
     * }
     * @endcode
     * and
     * @code{.bib}
     * @article{Dellinger1998,
     *   author  = {Dellinger, J. and Vasicek, D. and Sondergeld, C.},
     *   title   = {Kelvin notation for stabilizing elastic-constant inversion},
     *   journal = {Revue de l'Institut Fran{\c{c}}ais du P{\'e}trole},
     *   year    = {1998},
     *   volume  = {53},
     *   number  = {5},
     *   pages   = {709--719},
     *   url     = {http://sepwww.stanford.edu/oldsep/joe/Reprints/8IWSA.pdf},
     * }
     * @endcode
     * as well as the online reference found on
     * <a
     * href="https://en.wikipedia.org/wiki/Voigt_notation#Mandel_notation">this
     * wikipedia page</a> and <a
     * href="https://github.com/dealii/dealii/tree/master/tests/physics/notation-kelvin_02.cc">the
     * unit tests</a>.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    namespace Kelvin
    {
      /**
       * Input matrix has incorrect number of rows.
       */
      DeclException3(ExcNotationExcFullMatrixToTensorRowSize2,
                     int,
                     int,
                     int,
                     << "The number of rows in the input matrix is " << arg1
                     << ", but needs to be either " << arg2 << " or " << arg3
                     << ".");


      /**
       * Input matrix has incorrect number of rows.
       */
      DeclException4(ExcNotationExcFullMatrixToTensorRowSize3,
                     int,
                     int,
                     int,
                     int,
                     << "The number of rows in the input matrix is " << arg1
                     << ", but needs to be either " << arg2 << "," << arg3
                     << ", or " << arg4 << ".");


      /**
       * Input matrix has incorrect number of columns.
       */
      DeclException3(ExcNotationExcFullMatrixToTensorColSize2,
                     int,
                     int,
                     int,
                     << "The number of columns in the input matrix is " << arg1
                     << ", but needs to be either " << arg2 << " or " << arg3
                     << ".");


      /**
       * Input matrix has incorrect number of columns.
       */
      DeclException4(ExcNotationExcFullMatrixToTensorColSize3,
                     int,
                     int,
                     int,
                     int,
                     << "The number of columns in the input matrix is " << arg1
                     << ", but needs to be either " << arg2 << "," << arg3
                     << ", or " << arg4 << ".");


      /**
       * @name Forward operation: Tensor notation to Kelvin notation
       */
      //@{

      /**
       * Convert a scalar value to its compressed vector equivalent.
       *
       * The output vector has one entry.
       */
      template <typename Number>
      Vector<Number>
      to_vector(const Number &s);


      /**
       * Convert a rank-0 tensor to its compressed vector equivalent.
       *
       * The output vector has one entry.
       */
      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<0, dim, Number> &s);


      /**
       * Convert a rank-1 tensor to its compressed vector equivalent.
       *
       * The output vector has $dim$ entries.
       */
      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<1, dim, Number> &v);


      /**
       * Convert a rank-2 tensor to its compressed vector equivalent.
       *
       * The output vector has Tensor<2,dim>::n_independent_components entries.
       */
      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<2, dim, Number> &t);


      /**
       * Convert a rank-2 symmetric tensor to its compressed vector equivalent.
       *
       * The output vector has SymmetricTensor<2,dim>::n_independent_components
       * entries.
       */
      template <int dim, typename Number>
      Vector<Number>
      to_vector(const SymmetricTensor<2, dim, Number> &st);


      /**
       * Convert a scalar value to its compressed matrix equivalent.
       *
       * The output matrix will have one row and one column.
       */
      template <typename Number>
      FullMatrix<Number>
      to_matrix(const Number &s);


      /**
       * Convert a rank-0 tensor to its compressed matrix equivalent.
       *
       * The output matrix will have one row and one column.
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<0, dim, Number> &s);


      /**
       * Convert a rank-1 tensor to its compressed matrix equivalent.
       *
       * The output matrix will have $dim$ rows and one column.
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<1, dim, Number> &v);


      /**
       * Convert a rank-2 tensor to its compressed matrix equivalent.
       *
       * The output matrix will have $dim$ rows and $dim$ columns.
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<2, dim, Number> &t);


      /**
       * Convert a rank-2 symmetric tensor to its compressed matrix equivalent.
       *
       * The output matrix will have $dim$ rows and $dim$ columns, with the same
       * format as the equivalent function for non-symmetric tensors. This is
       * because it is not possible to compress the
       * SymmetricTensor<2,dim>::n_independent_components unique entries into a
       * square matrix.
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const SymmetricTensor<2, dim, Number> &st);

      /**
       * Convert a rank-3 tensor to its compressed matrix equivalent.
       *
       * The template arguments @p SubTensor1 and @p SubTensor2 determine how
       * the unrolling occurs, in particular how the elements of the rank-3
       * tensor are to be interpreted.
       *
       * So, for example, with the following two conversions
       * @code
       * Tensor<3,dim> r3_tnsr;      // All elements filled differently
       * Tensor<3,dim> r3_symm_tnsr; // Some elements filled symmetrically
       *
       * const FullMatrix<double> mtrx_1 =
       *   Physics::Notation::to_matrix<dim,
       *                                Tensor<2,dim>,
       *                                Tensor<1,dim>*>(r3_tnsr);
       * const FullMatrix<double> mtrx_2 =
       *   Physics::Notation::to_matrix<dim,
       *                                Tensor<1,dim>,
       *                                SymmetricTensor<2,dim>*>(r3_symm_tnsr);
       * @endcode
       * the matrix @p mtrx_1 will have $dim \times dim$ rows and $dim$ columns
       * (i.e. size Tensor<2,dim>::n_independent_components $\times$
       * Tensor<1,dim>::n_independent_components),
       * while those of the matrix @p mtrx_2 will have $dim$ rows and
       * $(dim \times dim + dim)/2$ columns
       * (i.e. size Tensor<1,dim>::n_independent_components $\times$
       * SymmetricTensor<2,dim>::n_independent_components), as it is assumed
       * that the entries corresponding to the alternation of the second and
       * third indices are equal. That is to say that
       * <code>r3_symm_tnsr[i][j][k] == r3_symm_tnsr[i][k][j]</code>.
       */
      template <int dim,
                typename SubTensor1 = Tensor<2, dim>,
                typename SubTensor2 = Tensor<1, dim>,
                typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<3, dim, Number> &t);


      /**
       * Convert a rank-4 tensor to its compressed matrix equivalent.
       *
       * The output matrix will have Tensor<2,dim>::n_independent_components
       * rows and Tensor<2,dim>::n_independent_components columns.
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<4, dim, Number> &t);


      /**
       * Convert a rank-4 symmetric tensor to its compressed matrix equivalent.
       *
       * The output matrix will have
       * SymmetricTensor<2,dim>::n_independent_components rows and
       * SymmetricTensor<2,dim>::n_independent_components columns.
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const SymmetricTensor<4, dim, Number> &st);

      //@}

      /**
       * @name Reverse operation: Kelvin notation to tensor notation
       */
      //@{

      /**
       * Convert a compressed vector to its equivalent scalar value.
       */
      template <typename Number>
      void
      to_tensor(const Vector<Number> &vec, Number &s);


      /**
       * Convert a compressed vector to its equivalent rank-0 tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<0, dim, Number> &s);


      /**
       * Convert a compressed vector to its equivalent rank-1 tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<1, dim, Number> &v);


      /**
       * Convert a compressed vector to its equivalent rank-2 tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<2, dim, Number> &t);


      /**
       * Convert a compressed vector to its equivalent rank-2 symmetric tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, SymmetricTensor<2, dim, Number> &st);


      /**
       * Convert a compressed matrix to its equivalent scalar value.
       */
      template <typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Number &s);


      /**
       * Convert a compressed matrix to its equivalent rank-0 tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<0, dim, Number> &s);


      /**
       * Convert a compressed matrix to its equivalent rank-1 tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<1, dim, Number> &v);


      /**
       * Convert a compressed matrix to its equivalent rank-2 tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<2, dim, Number> &t);


      /**
       * Convert a compressed matrix to its equivalent rank-2 symmetric tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &       mtrx,
                SymmetricTensor<2, dim, Number> &st);


      /**
       * Convert a compressed matrix to its equivalent rank-3 tensor.
       *
       * @note Based on the size of the matrix @p mtrx, some of the
       * components of @p t may be interpreted as having symmetric
       * counterparts. This is the reverse of the operation explained
       * in the documentation of the counterpart to_matrix()
       * function.
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<3, dim, Number> &t);


      /**
       * Convert a compressed matrix to its equivalent rank-4 tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<4, dim, Number> &t);


      /**
       * Convert a compressed matrix to its equivalent rank-4 symmetric tensor.
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &       mtrx,
                SymmetricTensor<4, dim, Number> &st);


      /**
       * A generic helper function that will convert a compressed vector
       * to its equivalent @p TensorType.
       */
      template <typename TensorType, typename Number>
      TensorType
      to_tensor(const Vector<Number> &vec);


      /**
       * A generic helper function that will convert a compressed matrix
       * to its equivalent @p TensorType.
       */
      template <typename TensorType, typename Number>
      TensorType
      to_tensor(const FullMatrix<Number> &vec);
      //@}

    } // namespace Kelvin

  } // namespace Notation
} // namespace Physics


#ifndef DOXYGEN


// ------------------------- inline functions ------------------------


namespace Physics
{
  namespace Notation
  {
    namespace Kelvin
    {
      namespace internal
      {
        /**
         * Return the tensor indices <code><row, column></code>
         * associated with a condensed component index. The
         * @p symmetric flag indicates whether or not the
         * @p component_n index is associated with a tensor that
         * has symmetric entries.
         */
        template <int dim>
        std::pair<unsigned int, unsigned int>
        indices_from_component(const unsigned int component_n,
                               const bool         symmetric);


        template <int dim>
        std::pair<unsigned int, unsigned int>
        indices_from_component(const unsigned int component_n, const bool)
        {
          AssertThrow(false, ExcNotImplemented());
          return std::make_pair(0u, 0u);
        }


        template <>
        inline std::pair<unsigned int, unsigned int>
        indices_from_component<1>(const unsigned int component_n, const bool)
        {
          AssertIndexRange(component_n, 1);

          return std::make_pair(0u, 0u);
        }


        template <>
        inline std::pair<unsigned int, unsigned int>
        indices_from_component<2>(const unsigned int component_n,
                                  const bool         symmetric)
        {
          if (symmetric == true)
            {
              Assert(
                (component_n < SymmetricTensor<2, 2>::n_independent_components),
                ExcIndexRange(component_n,
                              0,
                              SymmetricTensor<2, 2>::n_independent_components));
            }
          else
            {
              Assert((component_n < Tensor<2, 2>::n_independent_components),
                     ExcIndexRange(component_n,
                                   0,
                                   Tensor<2, 2>::n_independent_components));
            }

          static const unsigned int indices[4][2] = {{0, 0},
                                                     {1, 1},
                                                     {0, 1},
                                                     {1, 0}};
          return std::make_pair(indices[component_n][0],
                                indices[component_n][1]);
        }


        template <>
        inline std::pair<unsigned int, unsigned int>
        indices_from_component<3>(const unsigned int component_n,
                                  const bool         symmetric)
        {
          if (symmetric == true)
            {
              Assert(
                (component_n < SymmetricTensor<2, 3>::n_independent_components),
                ExcIndexRange(component_n,
                              0,
                              SymmetricTensor<2, 3>::n_independent_components));
            }
          else
            {
              Assert((component_n < Tensor<2, 3>::n_independent_components),
                     ExcIndexRange(component_n,
                                   0,
                                   Tensor<2, 3>::n_independent_components));
            }

          static const unsigned int indices[9][2] = {{0, 0},
                                                     {1, 1},
                                                     {2, 2},
                                                     {1, 2},
                                                     {0, 2},
                                                     {0, 1},
                                                     {1, 0},
                                                     {2, 0},
                                                     {2, 1}};
          return std::make_pair(indices[component_n][0],
                                indices[component_n][1]);
        }


        /**
         * Return the scaling factor to be applied to the
         * entry in the condensed vector.
         */
        template <int dim>
        double
        vector_component_factor(const unsigned int component_i,
                                const bool         symmetric)
        {
          if (symmetric == false)
            return 1.0;

          if (component_i < dim)
            return 1.0;
          else
            return numbers::SQRT2;
        }


        /**
         * Return the scaling factor to be applied to the
         * entry in the condensed matrix.
         */
        template <int dim>
        double
        matrix_component_factor(const unsigned int component_i,
                                const unsigned int component_j,
                                const bool         symmetric)
        {
          if (symmetric == false)
            return 1.0;

          // This case check returns equivalent of this result:
          // internal::vector_component_factor<dim>(component_i,symmetric)*internal::vector_component_factor<dim>(component_j,symmetric);
          if (component_i < dim && component_j < dim)
            return 1.0;
          else if (component_i >= dim && component_j >= dim)
            return 2.0;
          else // ((component_i >= dim && component_j < dim) || (component_i <
               // dim && component_j >= dim))
            return numbers::SQRT2;
        }

      } // namespace internal


      template <typename Number>
      Vector<Number>
      to_vector(const Number &s)
      {
        Vector<Number>     out(1);
        const unsigned int n_rows = out.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          out(r) = s;
        return out;
      }


      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<0, dim, Number> &s)
      {
        return to_vector(s.operator const Number &());
      }


      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<1, dim, Number> &v)
      {
        Vector<Number>     out(v.n_independent_components);
        const unsigned int n_rows = out.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            const unsigned int i = indices.first;
            out(r)               = v[i];
          }
        return out;
      }


      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<2, dim, Number> &t)
      {
        Vector<Number>     out(t.n_independent_components);
        const unsigned int n_rows = out.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second < dim, ExcInternalError());
            const unsigned int i = indices.first;
            const unsigned int j = indices.second;
            out(r)               = t[i][j];
          }
        return out;
      }


      template <int dim, typename Number>
      Vector<Number>
      to_vector(const SymmetricTensor<2, dim, Number> &st)
      {
        Vector<Number>     out(st.n_independent_components);
        const unsigned int n_rows = out.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, true);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second < dim, ExcInternalError());
            Assert(indices.second >= indices.first, ExcInternalError());
            const unsigned int i = indices.first;
            const unsigned int j = indices.second;

            const double factor =
              internal::vector_component_factor<dim>(r, true);

            out(r) = factor * st[i][j];
          }
        return out;
      }


      template <typename Number>
      FullMatrix<Number>
      to_matrix(const Number &s)
      {
        FullMatrix<Number> out(1, 1);
        out(0, 0) = s;
        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<0, dim, Number> &s)
      {
        return to_matrix(s.operator const Number &());
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<1, dim, Number> &v)
      {
        FullMatrix<Number> out(v.n_independent_components, 1);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            const unsigned int i = indices.first;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                Assert(c < 1, ExcInternalError());
                out(r, c) = v[i];
              }
          }
        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<2, dim, Number> &t)
      {
        FullMatrix<Number> out(dim, dim);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_i =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_i.first < dim, ExcInternalError());
            Assert(indices_i.second < dim, ExcInternalError());
            const unsigned int i = indices_i.first;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_j =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_j.first < dim, ExcInternalError());
                Assert(indices_j.second < dim, ExcInternalError());
                const unsigned int j = indices_j.second;

                out(r, c) = t[i][j];
              }
          }
        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const SymmetricTensor<2, dim, Number> &st)
      {
        return to_matrix(Tensor<2, dim, Number>(st));
      }


      namespace internal
      {
        template <typename TensorType>
        struct is_rank_2_symmetric_tensor : std::false_type
        {};

        template <int dim, typename Number>
        struct is_rank_2_symmetric_tensor<SymmetricTensor<2, dim, Number>>
          : std::true_type
        {};
      } // namespace internal


      template <int dim,
                typename SubTensor1,
                typename SubTensor2,
                typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<3, dim, Number> &t)
      {
        static_assert(
          (SubTensor1::dimension == dim && SubTensor2::dimension == dim),
          "Sub-tensor spatial dimension is different from those of the input tensor.");

        static_assert(
          (SubTensor1::rank == 2 && SubTensor2::rank == 1) ||
            (SubTensor1::rank == 1 && SubTensor2::rank == 2),
          "Cannot build a rank 3 tensor from the given combination of sub-tensors.");

        FullMatrix<Number> out(SubTensor1::n_independent_components,
                               SubTensor2::n_independent_components);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();

        if (SubTensor1::rank == 2 && SubTensor2::rank == 1)
          {
            const bool subtensor_is_rank_2_symmetric_tensor =
              internal::is_rank_2_symmetric_tensor<SubTensor1>::value;

            for (unsigned int r = 0; r < n_rows; ++r)
              {
                const std::pair<unsigned int, unsigned int> indices_ij =
                  internal::indices_from_component<dim>(
                    r, subtensor_is_rank_2_symmetric_tensor);
                Assert(indices_ij.first < dim, ExcInternalError());
                Assert(indices_ij.second < dim, ExcInternalError());
                if (subtensor_is_rank_2_symmetric_tensor)
                  {
                    Assert(indices_ij.second >= indices_ij.first,
                           ExcInternalError());
                  }
                const unsigned int i = indices_ij.first;
                const unsigned int j = indices_ij.second;

                const double factor = internal::vector_component_factor<dim>(
                  r, subtensor_is_rank_2_symmetric_tensor);

                for (unsigned int c = 0; c < n_cols; ++c)
                  {
                    const std::pair<unsigned int, unsigned int> indices_k =
                      internal::indices_from_component<dim>(c, false);
                    Assert(indices_k.first < dim, ExcInternalError());
                    const unsigned int k = indices_k.first;

                    if (subtensor_is_rank_2_symmetric_tensor)
                      out(r, c) = factor * t[i][j][k];
                    else
                      out(r, c) = t[i][j][k];
                  }
              }
          }
        else if (SubTensor1::rank == 1 && SubTensor2::rank == 2)
          {
            const bool subtensor_is_rank_2_symmetric_tensor =
              internal::is_rank_2_symmetric_tensor<SubTensor2>::value;

            for (unsigned int r = 0; r < n_rows; ++r)
              {
                const std::pair<unsigned int, unsigned int> indices_k =
                  internal::indices_from_component<dim>(r, false);
                Assert(indices_k.first < dim, ExcInternalError());
                const unsigned int k = indices_k.first;

                for (unsigned int c = 0; c < n_cols; ++c)
                  {
                    const std::pair<unsigned int, unsigned int> indices_ij =
                      internal::indices_from_component<dim>(
                        c, subtensor_is_rank_2_symmetric_tensor);
                    Assert(indices_ij.first < dim, ExcInternalError());
                    Assert(indices_ij.second < dim, ExcInternalError());
                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        Assert(indices_ij.second >= indices_ij.first,
                               ExcInternalError());
                      }
                    const unsigned int i = indices_ij.first;
                    const unsigned int j = indices_ij.second;

                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        const double factor =
                          internal::vector_component_factor<dim>(
                            c, subtensor_is_rank_2_symmetric_tensor);
                        out(r, c) = factor * t[k][i][j];
                      }
                    else
                      out(r, c) = t[k][i][j];
                  }
              }
          }
        else
          {
            AssertThrow(false, ExcNotImplemented());
          }

        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<4, dim, Number> &t)
      {
        FullMatrix<Number> out(
          Tensor<2, dim, Number>::n_independent_components,
          Tensor<2, dim, Number>::n_independent_components);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_ij =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_ij.first < dim, ExcInternalError());
            Assert(indices_ij.second < dim, ExcInternalError());
            const unsigned int i = indices_ij.first;
            const unsigned int j = indices_ij.second;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_kl =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_kl.first < dim, ExcInternalError());
                Assert(indices_kl.second < dim, ExcInternalError());
                const unsigned int k = indices_kl.first;
                const unsigned int l = indices_kl.second;

                out(r, c) = t[i][j][k][l];
              }
          }
        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const SymmetricTensor<4, dim, Number> &st)
      {
        FullMatrix<Number> out(
          SymmetricTensor<2, dim, Number>::n_independent_components,
          SymmetricTensor<2, dim, Number>::n_independent_components);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_ij =
              internal::indices_from_component<dim>(r, true);
            Assert(indices_ij.first < dim, ExcInternalError());
            Assert(indices_ij.second < dim, ExcInternalError());
            Assert(indices_ij.second >= indices_ij.first, ExcInternalError());
            const unsigned int i = indices_ij.first;
            const unsigned int j = indices_ij.second;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_kl =
                  internal::indices_from_component<dim>(c, true);
                Assert(indices_kl.first < dim, ExcInternalError());
                Assert(indices_kl.second < dim, ExcInternalError());
                Assert(indices_kl.second >= indices_kl.first,
                       ExcInternalError());
                const unsigned int k = indices_kl.first;
                const unsigned int l = indices_kl.second;

                const double factor =
                  internal::matrix_component_factor<dim>(r, c, true);

                out(r, c) = factor * st[i][j][k][l];
              }
          }
        return out;
      }


      template <typename Number>
      void
      to_tensor(const Vector<Number> &vec, Number &s)
      {
        Assert(vec.size() == 1, ExcDimensionMismatch(vec.size(), 1));
        s = vec(0);
      }


      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<0, dim, Number> &s)
      {
        return to_tensor(vec, s.operator Number &());
      }


      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<1, dim, Number> &v)
      {
        Assert(vec.size() == v.n_independent_components,
               ExcDimensionMismatch(vec.size(), v.n_independent_components));
        const unsigned int n_rows = vec.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            const unsigned int i = indices.first;
            v[i]                 = vec(r);
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<2, dim, Number> &t)
      {
        Assert(vec.size() == t.n_independent_components,
               ExcDimensionMismatch(vec.size(), t.n_independent_components));
        const unsigned int n_rows = vec.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second < dim, ExcInternalError());
            const unsigned int i = indices.first;
            const unsigned int j = indices.second;
            t[i][j]              = vec(r);
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, SymmetricTensor<2, dim, Number> &st)
      {
        Assert(vec.size() == st.n_independent_components,
               ExcDimensionMismatch(vec.size(), st.n_independent_components));
        const unsigned int n_rows = vec.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, true);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second < dim, ExcInternalError());
            Assert(indices.second >= indices.first, ExcInternalError());
            const unsigned int i = indices.first;
            const unsigned int j = indices.second;

            const double inv_factor =
              1.0 / internal::vector_component_factor<dim>(r, true);

            st[i][j] = inv_factor * vec(r);
          }
      }


      template <typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Number &s)
      {
        Assert(mtrx.m() == 1, ExcDimensionMismatch(mtrx.m(), 1));
        Assert(mtrx.n() == 1, ExcDimensionMismatch(mtrx.n(), 1));
        Assert(mtrx.n_elements() == 1,
               ExcDimensionMismatch(mtrx.n_elements(), 1));
        s = mtrx(0, 0);
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<0, dim, Number> &s)
      {
        return to_tensor(mtrx, s.operator Number &());
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<1, dim, Number> &v)
      {
        Assert(mtrx.m() == dim, ExcDimensionMismatch(mtrx.m(), dim));
        Assert(mtrx.n() == 1, ExcDimensionMismatch(mtrx.n(), 1));
        Assert(mtrx.n_elements() == v.n_independent_components,
               ExcDimensionMismatch(mtrx.n_elements(),
                                    v.n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second == 0, ExcInternalError());
            const unsigned int i = indices.first;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                Assert(c < 1, ExcInternalError());
                v[i] = mtrx(r, c);
              }
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<2, dim, Number> &t)
      {
        Assert(mtrx.m() == dim, ExcDimensionMismatch(mtrx.m(), dim));
        Assert(mtrx.n() == dim, ExcDimensionMismatch(mtrx.n(), dim));
        Assert(mtrx.n_elements() == t.n_independent_components,
               ExcDimensionMismatch(mtrx.n_elements(),
                                    t.n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_i =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_i.first < dim, ExcInternalError());
            Assert(indices_i.second < dim, ExcInternalError());
            const unsigned int i = indices_i.first;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_j =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_j.first < dim, ExcInternalError());
                Assert(indices_j.second < dim, ExcInternalError());
                const unsigned int j = indices_j.second;

                t[i][j] = mtrx(r, c);
              }
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &       mtrx,
                SymmetricTensor<2, dim, Number> &st)
      {
        // Its impossible to fit the (dim^2 + dim)/2 entries into a square
        // matrix We therefore assume that its been converted to a standard
        // tensor format using to_matrix (SymmetricTensor<2,dim,Number>) at some
        // point...
        Assert(mtrx.m() == dim, ExcDimensionMismatch(mtrx.m(), dim));
        Assert(mtrx.n() == dim, ExcDimensionMismatch(mtrx.n(), dim));
        Assert((mtrx.n_elements() ==
                Tensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.n_elements(),
                 Tensor<2, dim, Number>::n_independent_components));

        Tensor<2, dim, Number> tmp;
        to_tensor(mtrx, tmp);
        st = symmetrize(tmp);
        Assert((Tensor<2, dim, Number>(st) - tmp).norm() < 1e-12,
               ExcMessage(
                 "The entries stored inside the matrix were not symmetric"));
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<3, dim, Number> &t)
      {
        Assert((mtrx.m() == Tensor<1, dim, Number>::n_independent_components) ||
                 (mtrx.m() ==
                  Tensor<2, dim, Number>::n_independent_components) ||
                 (mtrx.m() ==
                  SymmetricTensor<2, dim, Number>::n_independent_components),
               ExcNotationExcFullMatrixToTensorColSize3(
                 mtrx.m(),
                 Tensor<1, dim, Number>::n_independent_components,
                 Tensor<2, dim, Number>::n_independent_components,
                 SymmetricTensor<2, dim, Number>::n_independent_components));
        Assert((mtrx.n() == Tensor<1, dim, Number>::n_independent_components) ||
                 (mtrx.n() ==
                  Tensor<2, dim, Number>::n_independent_components) ||
                 (mtrx.n() ==
                  SymmetricTensor<2, dim, Number>::n_independent_components),
               ExcNotationExcFullMatrixToTensorColSize3(
                 mtrx.n(),
                 Tensor<1, dim, Number>::n_independent_components,
                 Tensor<2, dim, Number>::n_independent_components,
                 SymmetricTensor<2, dim, Number>::n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        if (mtrx.n() == Tensor<1, dim, Number>::n_independent_components)
          {
            Assert(
              (mtrx.m() == Tensor<2, dim, Number>::n_independent_components) ||
                (mtrx.m() ==
                 SymmetricTensor<2, dim, Number>::n_independent_components),
              ExcNotationExcFullMatrixToTensorRowSize2(
                mtrx.m(),
                Tensor<2, dim, Number>::n_independent_components,
                SymmetricTensor<2, dim, Number>::n_independent_components));

            const bool subtensor_is_rank_2_symmetric_tensor =
              (mtrx.m() ==
               SymmetricTensor<2, dim, Number>::n_independent_components);

            for (unsigned int r = 0; r < n_rows; ++r)
              {
                const std::pair<unsigned int, unsigned int> indices_ij =
                  internal::indices_from_component<dim>(
                    r, subtensor_is_rank_2_symmetric_tensor);
                Assert(indices_ij.first < dim, ExcInternalError());
                Assert(indices_ij.second < dim, ExcInternalError());
                if (subtensor_is_rank_2_symmetric_tensor)
                  {
                    Assert(indices_ij.second >= indices_ij.first,
                           ExcInternalError());
                  }
                const unsigned int i = indices_ij.first;
                const unsigned int j = indices_ij.second;

                const double inv_factor =
                  1.0 / internal::vector_component_factor<dim>(
                          r, subtensor_is_rank_2_symmetric_tensor);

                for (unsigned int c = 0; c < n_cols; ++c)
                  {
                    const std::pair<unsigned int, unsigned int> indices_k =
                      internal::indices_from_component<dim>(c, false);
                    Assert(indices_k.first < dim, ExcInternalError());
                    const unsigned int k = indices_k.first;

                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        t[i][j][k] = inv_factor * mtrx(r, c);
                        t[j][i][k] = t[i][j][k];
                      }
                    else
                      t[i][j][k] = mtrx(r, c);
                  }
              }
          }
        else
          {
            Assert(
              (mtrx.m() == Tensor<1, dim, Number>::n_independent_components),
              ExcDimensionMismatch(
                mtrx.m(), Tensor<1, dim, Number>::n_independent_components));
            Assert(
              (mtrx.n() == Tensor<2, dim, Number>::n_independent_components) ||
                (mtrx.n() ==
                 SymmetricTensor<2, dim, Number>::n_independent_components),
              ExcNotationExcFullMatrixToTensorColSize2(
                mtrx.n(),
                Tensor<2, dim, Number>::n_independent_components,
                SymmetricTensor<2, dim, Number>::n_independent_components));

            const bool subtensor_is_rank_2_symmetric_tensor =
              (mtrx.n() ==
               SymmetricTensor<2, dim, Number>::n_independent_components);

            for (unsigned int r = 0; r < n_rows; ++r)
              {
                const std::pair<unsigned int, unsigned int> indices_k =
                  internal::indices_from_component<dim>(r, false);
                Assert(indices_k.first < dim, ExcInternalError());
                const unsigned int k = indices_k.first;

                for (unsigned int c = 0; c < n_cols; ++c)
                  {
                    const std::pair<unsigned int, unsigned int> indices_ij =
                      internal::indices_from_component<dim>(
                        c, subtensor_is_rank_2_symmetric_tensor);
                    Assert(indices_ij.first < dim, ExcInternalError());
                    Assert(indices_ij.second < dim, ExcInternalError());
                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        Assert(indices_ij.second >= indices_ij.first,
                               ExcInternalError());
                      }
                    const unsigned int i = indices_ij.first;
                    const unsigned int j = indices_ij.second;

                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        const double inv_factor =
                          1.0 / internal::vector_component_factor<dim>(
                                  c, subtensor_is_rank_2_symmetric_tensor);
                        t[k][i][j] = inv_factor * mtrx(r, c);
                        t[k][j][i] = t[k][i][j];
                      }
                    else
                      t[k][i][j] = mtrx(r, c);
                  }
              }
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<4, dim, Number> &t)
      {
        Assert((mtrx.m() == Tensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.m(), Tensor<2, dim, Number>::n_independent_components));
        Assert((mtrx.n() == Tensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.n(), Tensor<2, dim, Number>::n_independent_components));
        Assert(mtrx.n_elements() == t.n_independent_components,
               ExcDimensionMismatch(mtrx.n_elements(),
                                    t.n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_ij =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_ij.first < dim, ExcInternalError());
            Assert(indices_ij.second < dim, ExcInternalError());
            const unsigned int i = indices_ij.first;
            const unsigned int j = indices_ij.second;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_kl =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_kl.first < dim, ExcInternalError());
                Assert(indices_kl.second < dim, ExcInternalError());
                const unsigned int k = indices_kl.first;
                const unsigned int l = indices_kl.second;

                t[i][j][k][l] = mtrx(r, c);
              }
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &       mtrx,
                SymmetricTensor<4, dim, Number> &st)
      {
        Assert((mtrx.m() ==
                SymmetricTensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.m(),
                 SymmetricTensor<2, dim, Number>::n_independent_components));
        Assert((mtrx.n() ==
                SymmetricTensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.n(),
                 SymmetricTensor<2, dim, Number>::n_independent_components));
        Assert(mtrx.n_elements() == st.n_independent_components,
               ExcDimensionMismatch(mtrx.n_elements(),
                                    st.n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_ij =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_ij.first < dim, ExcInternalError());
            Assert(indices_ij.second < dim, ExcInternalError());
            const unsigned int i = indices_ij.first;
            const unsigned int j = indices_ij.second;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_kl =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_kl.first < dim, ExcInternalError());
                Assert(indices_kl.second < dim, ExcInternalError());
                const unsigned int k = indices_kl.first;
                const unsigned int l = indices_kl.second;

                const double inv_factor =
                  1.0 / internal::matrix_component_factor<dim>(r, c, true);

                st[i][j][k][l] = inv_factor * mtrx(r, c);
              }
          }
      }


      template <typename TensorType, typename Number>
      inline TensorType
      to_tensor(const Vector<Number> &vec)
      {
        TensorType out;
        to_tensor(vec, out);
        return out;
      }


      template <typename TensorType, typename Number>
      inline TensorType
      to_tensor(const FullMatrix<Number> &mtrx)
      {
        TensorType out;
        to_tensor(mtrx, out);
        return out;
      }

    } // namespace Kelvin
  }   // namespace Notation
} // namespace Physics


#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif
