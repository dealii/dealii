// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_vector_tools_integrate_difference_h
#define dealii_vector_tools_integrate_difference_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/vector_tools_common.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;
template <int dim>
class Quadrature;
template <int dim, int spacedim>
class Triangulation;
namespace hp
{
  template <int dim, int spacedim>
  class DoFHandler;
  template <int dim, int spacedim>
  class MappingCollection;
  template <int dim>
  class QCollection;
} // namespace hp


namespace VectorTools
{
  /**
   * @name Evaluation of functions and errors
   */
  //@{

  /**
   * Compute the cellwise error of the finite element solution.  Integrate the
   * difference between a reference function which is given as a continuous
   * function object, and a finite element function. The result of this
   * function is the vector @p difference that contains one value per active
   * cell $K$ of the triangulation. Each of the values of this vector $d$
   * equals
   * @f{align*}{
   * d_K = \| u-u_h \|_X
   * @f}
   * where $X$ denotes the norm chosen and $u$ represents the exact solution.
   *
   * It is assumed that the number of components of the function @p
   * exact_solution matches that of the finite element used by @p dof.
   *
   * To compute a global error norm of a finite element solution, use
   * VectorTools::compute_global_error() with the output vector computed with
   * this function.
   *
   * @param[in] mapping The mapping that is used when integrating the
   * difference $u-u_h$.
   * @param[in] dof The DoFHandler object that describes the finite element
   * space in which the solution vector lives.
   * @param[in] fe_function A vector with nodal values representing the
   * numerical approximation $u_h$. This vector needs to correspond to the
   * finite element space represented by @p dof.
   * @param[in] exact_solution The exact solution that is used to compute the
   * error.
   * @param[out] difference The vector of values $d_K$ computed as above.
   * @param[in] q The quadrature formula used to approximate the integral
   * shown above. Note that some quadrature formulas are more useful than
   * other in integrating $u-u_h$. For example, it is known that the $Q_1$
   * approximation $u_h$ to the exact solution $u$ of a Laplace equation is
   * particularly accurate (in fact, superconvergent, i.e. accurate to higher
   * order) at the 4 Gauss points of a cell in 2d (or 8 points in 3d) that
   * correspond to a QGauss(2) object. Consequently, because a QGauss(2)
   * formula only evaluates the two solutions at these particular points,
   * choosing this quadrature formula may indicate an error far smaller than
   * it actually is.
   * @param[in] norm The norm $X$ shown above that should be computed. If the
   * norm is NormType::Hdiv_seminorm, then the finite element on which this
   * function is called needs to have at least dim vector components, and the
   * divergence will be computed on the first div components. This works, for
   * example, on the finite elements used for the mixed Laplace (step-20) and
   * the Stokes equations (step-22).
   * @param[in] weight The additional argument @p weight allows to evaluate
   * weighted norms.  The weight function may be scalar, establishing a
   * spatially variable weight in the domain for all components equally. This
   * may be used, for instance, to only integrate over parts of the domain.
   * The weight function may also be vector-valued, with as many components as
   * the finite element: Then, different components get different weights. A
   * typical application is when the error with respect to only one or a
   * subset of the solution variables is to be computed, in which case the
   * other components would have weight values equal to zero. The
   * ComponentSelectFunction class is particularly useful for this purpose as
   * it provides such a "mask" weight. The weight function is expected to be
   * positive, but negative values are not filtered. The default value of this
   * function, a null pointer, is interpreted as "no weighting function",
   * i.e., weight=1 in the whole domain for all vector components uniformly.
   * @param[in] exponent This value denotes the $p$ used in computing
   * $L^p$-norms and $W^{1,p}$-norms. The value is ignored if a @p norm other
   * than NormType::Lp_norm, NormType::W1p_norm, or NormType::W1p_seminorm
   * is chosen.
   *
   *
   * See the general documentation of this namespace for more information.
   *
   * @note If the integration here happens over the cells of a
   * parallel::distribute::Triangulation object, then this function computes
   * the vector elements $d_K$ for an output vector with as many cells as
   * there are active cells of the triangulation object of the current
   * processor. However, not all active cells are in fact locally owned: some
   * may be ghost or artificial cells (see
   * @ref GlossGhostCell "here"
   * and
   * @ref GlossArtificialCell "here").
   * The vector computed will, in the case of a distributed triangulation,
   * contain zeros for cells that are not locally owned. As a consequence, in
   * order to compute the <i>global</i> $L_2$ error (for example), the errors
   * from different processors need to be combined, see
   * VectorTools::compute_global_error().
   *
   * Instantiations for this template are provided for some vector types (see
   * the general documentation of the namespace), but only for InVectors as in
   * the documentation of the namespace, OutVector only Vector<double> and
   * Vector<float>.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference(
    const Mapping<dim, spacedim> &                           mapping,
    const DoFHandler<dim, spacedim> &                        dof,
    const InVector &                                         fe_function,
    const Function<spacedim, typename InVector::value_type> &exact_solution,
    OutVector &                                              difference,
    const Quadrature<dim> &                                  q,
    const NormType &                                         norm,
    const Function<spacedim, double> *                       weight   = nullptr,
    const double                                             exponent = 2.);

  /**
   * Call the integrate_difference() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim@>(1)</tt>.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference(
    const DoFHandler<dim, spacedim> &                        dof,
    const InVector &                                         fe_function,
    const Function<spacedim, typename InVector::value_type> &exact_solution,
    OutVector &                                              difference,
    const Quadrature<dim> &                                  q,
    const NormType &                                         norm,
    const Function<spacedim, double> *                       weight   = nullptr,
    const double                                             exponent = 2.);

  /**
   * Same as above for hp.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference(
    const hp::MappingCollection<dim, spacedim> &             mapping,
    const hp::DoFHandler<dim, spacedim> &                    dof,
    const InVector &                                         fe_function,
    const Function<spacedim, typename InVector::value_type> &exact_solution,
    OutVector &                                              difference,
    const hp::QCollection<dim> &                             q,
    const NormType &                                         norm,
    const Function<spacedim, double> *                       weight   = nullptr,
    const double                                             exponent = 2.);

  /**
   * Call the integrate_difference() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim@>(1)</tt>.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference(
    const hp::DoFHandler<dim, spacedim> &                    dof,
    const InVector &                                         fe_function,
    const Function<spacedim, typename InVector::value_type> &exact_solution,
    OutVector &                                              difference,
    const hp::QCollection<dim> &                             q,
    const NormType &                                         norm,
    const Function<spacedim, double> *                       weight   = nullptr,
    const double                                             exponent = 2.);

  /**
   * Compute the cellwise error of the finite element solution.  Integrate the
   * difference between a reference function which is given as a continuous
   * function object, and a finite element function. The result of this
   * function is the vector @p difference that contains one value per active
   * cell $K$ of the triangulation. Each of the values of this vector $d$
   * equals
   * @f{align*}{
   * d_K = \| u-u_h \|_X
   * @f}
   * where $X$ denotes the norm chosen and $u$ represents the exact solution.
   *
   * @deprecated Use integrate_difference(const Mapping<dim, spacedim> &, const DoFHandler<dim, spacedim> &, const InVector &, const Function<spacedim, typename InVector::value_type> &, OutVector &, const Quadrature<dim> &, const NormType &, const Function<spacedim, double> *, const double) instead.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  DEAL_II_DEPRECATED typename std::enable_if<
    !std::is_same<typename InVector::value_type, double>::value>::type
  integrate_difference(const Mapping<dim, spacedim> &    mapping,
                       const DoFHandler<dim, spacedim> & dof,
                       const InVector &                  fe_function,
                       const Function<spacedim, double> &exact_solution,
                       OutVector &                       difference,
                       const Quadrature<dim> &           q,
                       const NormType &                  norm,
                       const Function<spacedim, double> *weight   = nullptr,
                       const double                      exponent = 2.);

  /**
   * Call the integrate_difference() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim@>(1)</tt>.
   *
   * @deprecated Use integrate_difference(const DoFHandler<dim, spacedim> &, const InVector &, const Function<spacedim, typename InVector::value_type> &exact_solution, OutVector &, const Quadrature<dim> &, const NormType &, const Function<spacedim, double> *, const double) instead.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  DEAL_II_DEPRECATED typename std::enable_if<
    !std::is_same<typename InVector::value_type, double>::value>::type
  integrate_difference(const DoFHandler<dim, spacedim> & dof,
                       const InVector &                  fe_function,
                       const Function<spacedim, double> &exact_solution,
                       OutVector &                       difference,
                       const Quadrature<dim> &           q,
                       const NormType &                  norm,
                       const Function<spacedim, double> *weight   = nullptr,
                       const double                      exponent = 2.);

  /**
   * Same as above for hp.
   *
   * @deprecated Use integrate_difference(const hp::MappingCollection<dim, spacedim> &, const hp::DoFHandler<dim, spacedim> &, const InVector &, const Function<spacedim, typename InVector::value_type> &, OutVector &, const hp::QCollection<dim> &, const NormType &, const Function<spacedim, double> *, const double) instead.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  DEAL_II_DEPRECATED typename std::enable_if<
    !std::is_same<typename InVector::value_type, double>::value>::type
  integrate_difference(const hp::MappingCollection<dim, spacedim> &mapping,
                       const hp::DoFHandler<dim, spacedim> &       dof,
                       const InVector &                            fe_function,
                       const Function<spacedim, double> &exact_solution,
                       OutVector &                       difference,
                       const hp::QCollection<dim> &      q,
                       const NormType &                  norm,
                       const Function<spacedim, double> *weight   = nullptr,
                       const double                      exponent = 2.);

  /**
   * Call the integrate_difference() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim@>(1)</tt>.
   *
   * @deprecated Use integrate_difference(const hp::DoFHandler<dim, spacedim> &, const InVector &, const Function<spacedim, typename InVector::value_type> &, OutVector &, const hp::QCollection<dim> &, const NormType &, const Function<spacedim, double> *, const double) instead.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  DEAL_II_DEPRECATED typename std::enable_if<
    !std::is_same<typename InVector::value_type, double>::value>::type
  integrate_difference(const hp::DoFHandler<dim, spacedim> &dof,
                       const InVector &                     fe_function,
                       const Function<spacedim, double> &   exact_solution,
                       OutVector &                          difference,
                       const hp::QCollection<dim> &         q,
                       const NormType &                     norm,
                       const Function<spacedim, double> *   weight   = nullptr,
                       const double                         exponent = 2.);

  /**
   * Take a Vector @p cellwise_error of errors on each cell with
   * <tt>tria.n_active_cells()</tt> entries and return the global
   * error as given by @p norm.
   *
   * The @p cellwise_error vector is typically an output produced by
   * VectorTools::integrate_difference() and you normally want to supply the
   * same value for @p norm as you used in VectorTools::integrate_difference().
   *
   * If the given Triangulation is a parallel::TriangulationBase, entries
   * in @p cellwise_error that do not correspond to locally owned cells are
   * assumed to be 0.0 and a parallel reduction using MPI is done to compute
   * the global error.
   *
   * @param tria The Triangulation with active cells corresponding with the
   * entries in @p cellwise_error.
   * @param cellwise_error Vector of errors on each active cell.
   * @param norm The type of norm to compute.
   * @param exponent The exponent $p$ to use for $L^p$-norms and
   * $W^{1,p}$-norms. The value is ignored if a @p norm other
   * than NormType::Lp_norm, NormType::W1p_norm, or NormType::W1p_seminorm
   * is chosen.
   *
   * @note Instantiated for type Vector<double> and Vector<float>.
   */
  template <int dim, int spacedim, class InVector>
  double
  compute_global_error(const Triangulation<dim, spacedim> &tria,
                       const InVector &                    cellwise_error,
                       const NormType &                    norm,
                       const double                        exponent = 2.);

  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_integrate_difference_h
