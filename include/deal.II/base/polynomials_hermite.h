// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#ifndef dealii_polynomials_hermite_h
#define dealii_polynomials_hermite_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>

#include <vector>



DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Polynomials
 * @{
 */

namespace Polynomials
{
  /**
   * This class implements Hermite interpolation polynomials (see
   * @cite CiarletRiavart1972interpolation) enforcing the maximum
   * possible level of regularity $r$ in the FEM basis given a
   * polynomial degree of $2r+1$. The polynomials all represent
   * either a non-zero shape value or derivative at $x=0$ and $x=1$
   * on the reference interval $x \in [0,1]$.
   *
   * Indices $j = 0, 1, \dots, r$ refer to polynomials corresponding
   * to a non-zero derivative (or shape value for $j=0$) of
   * order $j$ at $x=0$, and indices $j = r+1, r+2, \dots, 2r+1$
   * refer to polynomials with a non-zero derivative of order
   * $j-(r+1)$ (or value for $j=r+1$) at $x=1$. In particular, the
   * $0^{th}$ function has a value of $1$ at $x=0$, and the
   * $(r+1)^{th}$ function has a value of $1$ at $x=1$.The basis is
   * rescaled such that a function corresponding to a non-zero $j^{th}$
   * derivative has derivative value $j! 4^{j}$ at the corresponding
   * node. This is done to prevent the $L^{2}$-norm of the basis functions
   * from reducing exponentially with the chosen regularity.
   */
  class PolynomialsHermite : public Polynomial<double>
  {
  public:
    /**
     * Constructor for an individual Hermite polynomial. We write $f_{j}$
     * for a polynomial that has a non-zero $j^{th}$ derivative at $x=0$
     * and $g_{j}$ for a polynomial with a non-zero $j^{th}$ derivative
     * at $x=1$, meaning $f_{j}$ will have @p index $=j$ and $g_{j}$ will
     * have @p index $= j + \mathtt{regularity} + 1$. The resulting
     * polynomials will be degree $2\times \mathtt{regularity} +1$
     * and obey the following conditions:
     * @f{align*}{
     * &\begin{matrix}
     *   \left. \frac{d^{i}}{dx^{i}} f_{j}(x) \right\vert_{x=0}
     *          = i! 4^{i} \delta_{i, j}, \hfill
     *          &\qquad \hfill 0 \leq i \leq \mathtt{regularity}, \\
     *   \left. \frac{d^{i}}{dx^{i}} f_{j}(x) \right\vert_{x=1}
     *          = 0, \hfill &\qquad \hfill 0 \leq i \leq \mathtt{regularity},
     * \end{matrix} \qquad 0 \leq j \leq \mathtt{regularity}, \\
     * &\begin{matrix}
     *  \left. \frac{d^{i}}{dx^{i}} g_{j}(x) \right\vert_{x=0}
     *          = 0, \hfill &\qquad \hfill 0 \leq i \leq \mathtt{regularity}, \\
     * \left. \frac{d^{i}}{dx^{i}} g_{j}(x) \right\vert_{x=1}
     *          = i! 4^{i} \delta_{i, j}, \hfill
     *          &\qquad \hfill 0 \leq i \leq \mathtt{regularity},
     * \end{matrix} \qquad 0 \leq j \leq \mathtt{regularity},
     * @f}
     * where $\delta_{i,j}$ is equal to $1$ whenever $i=j$,
     * and equal to $0$ otherwise. These polynomials have explicit
     * formulas given by
     * @f{align*}{
     *   f_{j}(x) &= 4^{j} x^{j} (1-x)^{\mathtt{regularity}+1}
     * \sum_{k=0}^{\mathtt{regularity} - j} \;^{\mathtt{regularity} + k} C_{k}
     * x^{k}, \\ g_{j}(x) &= 4^{j} x^{\mathtt{regularity}+1} (x-1)^{j}
     * \sum_{k=0}^{\mathtt{regularity} - j} \;^{\mathtt{regularity} + k} C_{k}
     * (1-x)^{k},
     * @f}
     * where $^{n} C_{r} = \frac{n!}{r!(n-r)!}$ is the $r^{th}$ binomial
     * coefficient of degree $n, \; 0 \leq r \leq n$.
     *
     * @param regularity The highest derivative for which the basis
     * is used to enforce regularity.
     * @param index The local index of the generated polynomial in the
     * Hermite basis.
     */
    PolynomialsHermite(const unsigned int regularity, const unsigned int index);

    /**
     * This function generates a vector of Polynomial objects
     * representing a complete basis of degree $2\times\mathtt{regularity} +1$
     * on the reference interval $[0,1]$.
     *
     * @param regularity The generated basis can be used to strongly
     * enforce continuity in all derivatives up to and including this
     * order.
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int regularity);

  protected:
    /**
     * Degree of the polynomial basis being used.
     */
    unsigned int degree;

    /**
     * The order of the highest derivative in which the Hermite
     * basis can be used to impose continuity across element
     * boundaries. It's related to the degree $p$ by
     * $p = 2 \times\mathtt{regularity} +1$.
     */
    unsigned int regularity;

    /**
     * This variable stores the derivative that the shape function
     * corresponds to at the element boundary given by <code>side</code>.
     */
    unsigned int side_index;

    /**
     * This stores whether the shape function corresponds to a non-zero
     * value or derivative at $x=0$ on the reference interval
     * ($\mathtt{side} =0$) or at $x=1$ ($\mathtt{side} =1$).
     */
    unsigned int side;
  };
} // namespace Polynomials
/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
