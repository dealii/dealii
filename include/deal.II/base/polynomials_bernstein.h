#include <deal.II/base/polynomial.h>
#include <fstream>
#include <iostream>

#ifndef dealii__polynomials_bernstein_h
#define dealii__polynomials_bernstein_h


DEAL_II_NAMESPACE_OPEN

/**
 * This class implements Bernstein basis polynomials of desire degree as
 * described in http://www.idav.ucdavis.edu/education/CAGDNotes/Bernstein-
 * Polynomials.pdf in the paragraph "Converting from the Bernstein Basis to
 * the Power Basis".
 *
 * They are used to create the Bernstein finite element FE_Bernstein.
 *
 * @ingroup Polynomials
 * @author Luca Heltai, Marco Tezzele
 * @date 2013, 2015
 */
template <typename number>
class PolynomialsBernstein : public Polynomials::Polynomial<number>
{
public:
  /**
   * Construct the @p index -th Bernstein Polynomial of degree @p degree.
   *
   * @f{align*}{
   * B_{\text{index}, \text{degree}} (t)
   *   &= \text{binom}(\text{degree}, \text{index})
   *      \cdot t^{\text{index}}
   *      \cdot (1 - t)^{\text{degree} - \text{index}} \\
   *   &= \sum_{i = \text{index}}^\text{degree}
   *      \cdot (-1)^{i - \text{index}}
   *      \cdot \text{binom}(\text{degree}, i)
   *      \cdot \text{binom}(i, \text{index})
   *      \cdot t^i
   * @f}
   *
   * @param index
   * @param degree
   */
  PolynomialsBernstein (
    const unsigned int index,
    const unsigned int degree);
};


template <typename number>
std::vector<Polynomials::Polynomial<number> >
generate_complete_bernstein_basis (
  const unsigned int degree)
{
  std::vector<Polynomials::Polynomial<number> > v;
  for (unsigned int i = 0; i < degree + 1; ++i)
    v.push_back(PolynomialsBernstein<number>(i, degree));
  return v;
}

DEAL_II_NAMESPACE_CLOSE

#endif
