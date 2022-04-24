// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#ifndef dealii_solver_gcr_h
#define dealii_solver_gcr_h

#include <deal.II/lac/solver_gmres.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Implementation of the GCR Method. The stopping criterion is the norm of the
 * residual.
 */
template <class VectorType = Vector<double>>
class SolverGCR : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    explicit AdditionalData(const unsigned int max_n_tmp_vectors    = 30,
                            const bool         use_default_residual = true);

    /**
     * Maximum number of temporary vectors. This parameter controls the size
     * of the Arnoldi basis, which for historical reasons is
     * #max_n_tmp_vectors-2. SolverGMRES assumes that there are at least three
     * temporary vectors, so this value must be greater than or equal to three.
     */
    unsigned int max_n_tmp_vectors;

    /**
     * Flag for the default residual that is used to measure convergence.
     */
    bool use_default_residual;
  };

  SolverGCR(SolverControl &       solver_control,
            const AdditionalData &data = AdditionalData());

  /**
   *  Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

private:
  /**
   * Maximum number of temporary vectors.
   */
  const unsigned int max_n_tmp_vectors;

  /**
   * Flag for the default residual that is used to measure convergence.
   */
  const bool use_default_residual;
};



template <class VectorType>
inline SolverGCR<VectorType>::AdditionalData::AdditionalData(
  const unsigned int max_n_tmp_vectors,
  const bool         use_default_residual)
  : max_n_tmp_vectors(max_n_tmp_vectors)
  , use_default_residual(use_default_residual)
{}



template <class VectorType>
inline SolverGCR<VectorType>::SolverGCR(SolverControl &       solver_control,
                                        const AdditionalData &data)
  : SolverBase<VectorType>(solver_control)
  , max_n_tmp_vectors(data.max_n_tmp_vectors)
  , use_default_residual(data.use_default_residual)
{
  solver_control.set_max_steps(max_n_tmp_vectors);
}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
inline void
SolverGCR<VectorType>::solve(const MatrixType &        A,
                             VectorType &              x,
                             const VectorType &        b,
                             const PreconditionerType &preconditioner)
{
  using number = typename VectorType::value_type;

  SolverControl::State conv = SolverControl::iterate;

  typename VectorMemory<VectorType>::Pointer search_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer Asearch_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer p_pointer(this->memory);

  VectorType &search  = *search_pointer;
  VectorType &Asearch = *Asearch_pointer;
  VectorType &p       = *p_pointer;

  std::vector<typename VectorType::value_type> Hn_preloc;
  Hn_preloc.reserve(max_n_tmp_vectors);

  internal::SolverGMRESImplementation::TmpVectors<VectorType> H_vec(
    max_n_tmp_vectors, this->memory);
  internal::SolverGMRESImplementation::TmpVectors<VectorType> Hd_vec(
    max_n_tmp_vectors, this->memory);

  search.reinit(x);
  Asearch.reinit(x);
  p.reinit(x);

  typename VectorType::value_type res = 0.0;

  A.vmult(p, x);
  p.add(-1., b);

  if (use_default_residual == true)
    res = p.l2_norm();

  preconditioner.vmult(search, p);

  if (use_default_residual == false)
    res = search.l2_norm();

  unsigned int it = 0;

  conv = this->iteration_status(it, res, x);
  if (conv != SolverControl::iterate)
    return;

  while (conv == SolverControl::iterate)
    {
      it++;

      H_vec(it - 1, x);
      Hd_vec(it - 1, x);

      Hn_preloc.resize(it);

      A.vmult(Asearch, search);

      for (unsigned int i = 0; i < it - 1; ++i)
        {
          const auto temptest = (H_vec[i] * Asearch) / Hn_preloc[i];
          Asearch.add(-temptest, H_vec[i]);
          search.add(-temptest, Hd_vec[i]);
        }

      const auto nAsearch_new = Asearch.norm_sqr();
      Hn_preloc[it - 1]       = nAsearch_new;
      H_vec[it - 1]           = Asearch;
      Hd_vec[it - 1]          = search;

      Assert(std::abs(nAsearch_new) != 0., ExcDivideByZero());

      const auto c_preloc = (Asearch * p) / nAsearch_new;
      x.add(-c_preloc, search);
      p.add(-c_preloc, Asearch);

      if (use_default_residual == true)
        res = p.l2_norm();

      preconditioner.vmult(search, p);

      if (use_default_residual == false)
        res = search.l2_norm();

      conv = this->iteration_status(it, res, x);
    }

  if (conv != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence(it, res));
}


DEAL_II_NAMESPACE_CLOSE

#endif
