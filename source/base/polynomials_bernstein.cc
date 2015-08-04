
#include <deal.II/base/polynomials_bernstein.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Binomial
{
  unsigned int
  factorial (
    unsigned int k)
  {
    if (k > 1)
      return k * factorial(k - 1);
    else
      return 1;
  }

  inline unsigned int
  binomial (
    const unsigned int nValue, const unsigned int nValue2)
  {
    Assert(nValue >= nValue2,
           ExcMessage("nValue should be greater or equal than nValue2"));
    if (nValue2 == 1)
      return nValue;
    else
      return ((factorial(nValue))
              / (factorial(nValue2) * factorial((nValue - nValue2))));
  }

  template <typename number>
  std::vector<number>
  get_bernstein_coefficients (
    const unsigned int k, const unsigned int n)
  {
    Assert(n>0, ExcMessage("Bernstein polynomial needs to be of degree > 0."));
    AssertIndexRange(k, n+1);
    std::vector<number> coeff(n + 1, number(0.0));
    for (unsigned int i = k; i < n + 1; ++i)
      coeff[i] = (pow(number(-1), number(i - k)) * binomial(n, i)
                  * binomial(i, k));
    return coeff;
  }

}

template <typename number>
PolynomialsBernstein<number>:: PolynomialsBernstein (
  const unsigned int index, const unsigned int degree)
  :
  Polynomials::Polynomial<number>(
   Binomial::get_bernstein_coefficients<number>(index, degree))
{
}


#include "polynomials_bernstein.inst"

DEAL_II_NAMESPACE_CLOSE
