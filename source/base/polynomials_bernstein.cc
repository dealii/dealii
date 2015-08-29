#include <deal.II/base/polynomials_bernstein.h>

#include <boost/math/special_functions/binomial.hpp>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace
{
  template <typename number>
  std::vector<number>
  get_bernstein_coefficients (
    const unsigned int k, const unsigned int n)
  {
    Assert(n>0, ExcMessage("Bernstein polynomial needs to be of degree > 0."));
    AssertIndexRange(k, n+1);
    std::vector<number> coeff(n + 1, number(0.0));
    for (unsigned int i = k; i < n + 1; ++i)
      {
        coeff[i] = ((i - k) % 2 == 0 ? 1 : -1)
                   * boost::math::binomial_coefficient<number>(n, i)
                   * boost::math::binomial_coefficient<number>(i, k);
      }
    return coeff;
  }
}

template <typename number>
PolynomialsBernstein<number>:: PolynomialsBernstein (
  const unsigned int index, const unsigned int degree)
  :
  Polynomials::Polynomial<number>(
    get_bernstein_coefficients<number>(index, degree))
{
}


#include "polynomials_bernstein.inst"

DEAL_II_NAMESPACE_CLOSE
