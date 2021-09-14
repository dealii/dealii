#pragma once

#include <deal.II/lac/solver.h>

namespace dealii
{
  template <typename VectorType>
  class SolverGCR : public SolverBase<VectorType>
  {
  public:
    SolverGCR(SolverControl &solver_control, const unsigned int GCRmaxit = 40)
      : SolverBase<VectorType>(solver_control)
      , GCRmaxit(GCRmaxit)
    {
      H_vec.reserve(GCRmaxit);
      Hd_vec.reserve(GCRmaxit);
    }

    template <typename MatrixType, typename PreconditionerType>
    void
    solve(const MatrixType &        A,
          VectorType &              x,
          const VectorType &        b,
          const PreconditionerType &preconditioner)
    {
      using number = typename VectorType::value_type;

      SolverControl::State conv = SolverControl::iterate;

      search.reinit(x);
      Asearch.reinit(x);
      p.reinit(x);

      A.vmult(p, x);
      p.add(-1., b);
      preconditioner.vmult(search, p);

      double res = search.l2_norm();

      unsigned int it = 0;

      conv = this->iteration_status(it, res, x);
      if (conv != SolverControl::iterate)
        return;

      while (conv == SolverControl::iterate)
        {
          AssertIndexRange(it, GCRmaxit);

          it++;

          if (H_vec.size() < it)
            {
              H_vec.resize(H_vec.size() + 1);
              Hd_vec.resize(Hd_vec.size() + 1);
              Hn_preloc.resize(Hn_preloc.size() + 1);

              H_vec.back().reinit(x);
              Hd_vec.back().reinit(x);
            }

          A.vmult(Asearch, search);

          for (unsigned int i = 0; i < it - 1; ++i)
            {
              const double temptest = (H_vec[i] * Asearch) / Hn_preloc[i];
              Asearch.add(-temptest, H_vec[i]);
              search.add(-temptest, Hd_vec[i]);
            }

          const double nAsearch_new = Asearch.norm_sqr();
          Hn_preloc[it - 1]         = nAsearch_new;
          H_vec[it - 1]             = Asearch;
          Hd_vec[it - 1]            = search;

          const double c_preloc = (Asearch * p) / nAsearch_new;
          x.add(c_preloc, search);
          p.add(-c_preloc, Asearch);

          preconditioner.vmult(search, p);
          res = search.l2_norm();

          conv = this->iteration_status(it, res, x);
        }

      if (conv != SolverControl::success)
        AssertThrow(false, SolverControl::NoConvergence(it, res));
    }

  private:
    const unsigned int GCRmaxit;

    mutable VectorType search;
    mutable VectorType Asearch;
    mutable VectorType p;

    mutable std::vector<typename VectorType::value_type> Hn_preloc;
    mutable std::vector<VectorType>                      H_vec;
    mutable std::vector<VectorType>                      Hd_vec;
  };

} // namespace dealii
