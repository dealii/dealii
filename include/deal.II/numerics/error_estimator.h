// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__error_estimator_h
#define __deal2__error_estimator_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/dofs/function_map.h>
#include <deal.II/fe/component_mask.h>
#include <map>

DEAL_II_NAMESPACE_OPEN


template <int, int> class DoFHandler;
template <int, int> class Mapping;
template <int> class Quadrature;

namespace hp
{
  template <int, int> class DoFHandler;
  template <int> class QCollection;
}



/**
 *  Implementation of the error indicator by Kelly, Gago, Zienkiewicz
 *  and Babuska.  This error indicator tries to approximate the error
 *  per cell by integration of the jump of the gradient of the
 *  solution along the faces of each cell.  It can be understood as a
 *  gradient recovery estimator; see the survey of Ainsworth for a
 *  complete discussion.
 *
 *  @note In spite of the name, this is not truly an a posteriori
 *  error estimator, even if applied to the Poisson problem only. It
 *  gives good hints for mesh refinement, but the estimate is not to
 *  be trusted. For higher order trial spaces the integrals computed
 *  here tend to zero faster than the error itself, thus ruling out
 *  the values as error estimators.
 *
 *  The error estimator really only estimates the error for the generalized
 *  Poisson equation $-\nabla\cdot a(x) \nabla u = f$ with either Dirichlet
 *  boundary conditions or generalized Neumann boundary conditions involving
 *  the conormal derivative $a\frac{du}{dn} = g$.
 *
 *  The error estimator returns a vector of estimated errors per cell which
 *  can be used to feed the <tt>Triangulation<dim>::refine_*</tt> functions. This
 *  vector contains elements of data type @p float, rather than @p double,
 *  since accuracy is not so important here, and since this can save rather
 *  a lot of memory, when using many cells.
 *
 *
 *  <h3>Implementation</h3>
 *
 *  In principle, the implementation of the error estimation is simple: let
 *  \f[
 *  \eta_K^2 = \frac h{24} \int_{\partial K} \left[a \frac{\partial u_h}{\partial n}\right]^2 do
 *  \f]
 *  be the error estimator for cell $K$. $[\cdot]$ denotes the jump of the
 *  argument at the face. In the paper of Ainsworth, $h$ is divided by $24$,
 *  but this factor is a bit esoteric, stemming from interpolation estimates
 *  and stability constants which may hold for the Poisson problem, but may
 *  not hold for more general situations. In the implementation, this factor
 *  is considered, but may lead to wrong results. You may scale the vector
 *  appropriately afterwards.
 *
 *  To perform the integration, use is made of the FEFaceValues and
 *  FESubfaceValues classes. The integration is performed by looping
 *  over all cells and integrating over faces that are not yet treated.
 *  This way we avoid integration on faces twice, once for each time we
 *  visit one of the adjacent cells. In a second loop over all cells, we
 *  sum up the contributions of the faces (which are the integrated
 *  square of the jumps) of each cell and take the square root.
 *
 *  The integration is done using a quadrature formula on the face.
 *  For linear trial functions (FEQ1), the QGauss2 or even the
 *  QMidpoint rule will suffice. For higher order elements, it is
 *  necessary to utilize higher order quadrature formulae as well.
 *
 *  We store the contribution of each face in a @p map, as provided by the
 *  C++ standard library, with the iterator pointing to that face being the
 *  key into the map. In fact, we do not store the indicator per face, but
 *  only the integral listed above. When looping the second time over all
 *  cells, we have to sum up the contributions of the faces, multiply them
 *  with $\frac h{24}$ and take the square root. By doing the multiplication
 *  with $h$ in the second loop, we avoid problems to decide with which $h$
 *  to multiply, that of the cell on the one or that of the cell on the other
 *  side of the face.
 *
 *  $h$ is taken to be the greatest length of the diagonals of the cell. For
 *  more or less uniform cells without deformed angles, this coincides with
 *  the diameter of the cell.
 *
 *
 *  <h3>Vector-valued functions</h3>
 *
 *  If the finite element field for which the error is to be estimated
 *  is vector-valued, i.e. the finite element has more than one
 *  component, then all the above can be applied to all or only some
 *  components at the same time. The main function of this class takes
 *  a list of flags denoting the components for which components the
 *  error estimator is to be applied; by default, it is a list of only
 *  @p trues, meaning that all variables shall be treated.
 *
 *  In case the different components of a field have different
 *  physical meaning (for example velocity and pressure in the Stokes
 *  equations), it would be senseless to use the same coefficient for
 *  all components. In that case, you can pass a function with as many
 *  components as there are components in the finite element field and
 *  each component of the error estimator will then be weighted by the
 *  respective component in this coefficient function. In the other
 *  case, when all components have the same meaning (for example the
 *  displacements in Lame's equations of elasticity), you can specify
 *  a scalar coefficient which will then be used for all components.
 *
 *
 *  <h3>%Boundary values</h3>
 *
 *  If the face is at the boundary, i.e. there is no neighboring cell to which
 *  the jump in the gradiend could be computed, there are two possibilities:
 *  <ul>
 *  <li> The face belongs to a Dirichlet boundary. Then the face is not
 *    considered, which can be justified looking at a dual problem technique and
 *    should hold exactly if the boundary can be approximated exactly by the
 *    finite element used (i.e. it is a linear boundary for linear finite elements,
 *    quadratic for isoparametric quadratic elements, etc). For boundaries which
 *    can not be exactly approximated, one should consider the difference
 *    $z-z_h$ on the face, $z$ being a dual problem's solution which is zero at
 *    the true boundary and $z_h$ being an approximation, which in most cases
 *    will be zero on the numerical boundary. Since on the numerical boundary
 *    $z$ will not be zero in general, we would get another term here, but this
 *    one is neglected for practical reasons, in the hope that the error made
 *    here will tend to zero faster than the energy error we wish to estimate.
 *
 *    Though no integration is necessary, in the list of face contributions we
 *    store a zero for this face, which makes summing up the contributions of
 *    the different faces to the cells easier.
 *
 *  <li> The face belongs to a Neumann boundary.  In this case, the
 *    contribution of the face $F\in\partial K$ looks like
 *    \f[ \int_F \left|g-a\frac{\partial u_h}{\partial n}\right|^2 ds \f]
 *    where $g$ is the Neumann boundary function. If the finite element is
 *    vector-valued, then obviously the function denoting the Neumann boundary
 *    conditions needs to be vector-valued as well.
 *
 *  <li> No other boundary conditions are considered.
 *  </ul>
 *
 *  In practice, if you have Robin boundary conditions or are too lazy to
 *  accurately describe Neumann values, then this is rarely an issue: if you
 *  don't say anything in the map about a particular part of the boundary then
 *  the Kelly indicator will simply assume that the solution is correct on
 *  that part of the boundary and not touch it. Of course, if you have a have
 *  a Neumann or Robin boundary, that isn't quite true, there is going to be a
 *  difference between the normal derivative of the numerical solution and the
 *  Neumann values these normal derivatives should equal. So if we simply
 *  ignore those parts of the boundary, we'll underestimate the error. In
 *  practice, this rarely appears to be a problem -- you may not refine the
 *  cell this time around but you'll probably refine it in the next refinement
 *  step and everything is good again. After all, for all problems but the
 *  Laplace equation, the Kelly indicator is only an indicator, not an
 *  estimator, and so the values it computes are not exact error
 *  representations anyway.
 *
 *
 *  <h3>Handling of hanging nodes</h3>
 *
 *  The integration along faces with hanging nodes is quite tricky, since one
 *  of the elements has to be shifted one level up or down. See the
 *  documentation for the FESubFaceValues class for more information about
 *  technical issues regarding this topic.
 *
 *  In praxi, since we integrate over each face only once, we do this when we
 *  are on the coarser one of the two cells adjacent to a subface (a subface
 *  is defined to be the child of a face; seen from the coarse cell, it is a
 *  subface, while seen from the refined cell it is one of its faces). The
 *  reason is that finding neighborship information is a bit easier then, but
 *  that's all practical reasoning, nothing fundamental.
 *
 *  Since we integrate from the coarse side of the face, we have the mother
 *  face readily at hand and store the result of the integration over that
 *  mother face (being the sum of the integrals along the subfaces) in the
 *  abovementioned map of integrals as well. This consumes some memory more
 *  than needed, but makes the summing up of the face contributions to the
 *  cells easier, since then we have the information from all faces of all
 *  cells at hand and need not think about explicitly determining whether
 *  a face was refined or not. The same applies for boundary faces, see
 *  above.
 *
 *
 *  <h3>Multiple solution vectors</h3>
 *
 *  In some cases, for example in time-dependent problems, one would
 *  like to compute the error estimates for several solution vectors
 *  on the same grid at once, with the same coefficients, boundary
 *  condition object, etc, e.g. for the solutions on several
 *  successive time steps. One could then call the functions of this
 *  class several times for each solution. However, the largest factor
 *  in the computation of the error estimates (in terms of computing
 *  time) is initialization of FEFaceValues and FESubFaceValues
 *  objects, and iterating through all faces and subfaces. If the
 *  solution vectors live on the same grid, this effort can be reduced
 *  significantly by treating all solution vectors at the same time,
 *  initializing the FEFaceValues objects only once per cell and for
 *  all solution vectors at once, and also only looping through the
 *  triangulation only once. For this reason, besides the @p estimate
 *  function in this class that takes a single input vector and
 *  returns a single output vector, there is also a function that
 *  accepts several in- and output vectors at the same time.
 *
 *  @ingroup numerics
 *  @author Wolfgang Bangerth, 1998, 1999, 2000, 2004, 2006; parallelization by Thomas Richter, 2000
 */
template <int dim, int spacedim=dim>
class KellyErrorEstimator
{
public:
  /**
   * Implementation of the error estimator described above. You may give a
   * coefficient, but there is a default value which denotes the constant
   * coefficient with value one. The coefficient function may either be a
   * scalar one, in which case it is used for all components of the finite
   * element, or a vector-valued one with as many components as there are in
   * the finite element; in the latter case, each component is weighted by the
   * respective component in the coefficient.
   *
   * You might give a list of components you want to evaluate, in case the
   * finite element used by the DoFHandler object is vector-valued. You then
   * have to set those entries to true in the bit-vector @p component_mask
   * (see @ref GlossComponentMask) for which the respective component is to be
   * used in the error estimator. The default is to use all components, which
   * is done by either providing a bit-vector with all-set entries, or an
   * empty bit-vector.
   *
   * The @p subdomain_id parameter indicates whether we shall compute
   * indicators for all cells (in case its value is the default,
   * <tt>numbers::invalid_unsigned_int</tt>), or only for the cells belonging
   * to a certain subdomain with the given indicator. The latter case is used
   * for parallel computations where all processor nodes have the global grid
   * stored, and could well compute all the indicators for all cells
   * themselves, but where it is more efficient to have each process compute
   * only indicators for the cells it owns, and have them exchange the
   * resulting information afterwards. This is in particular true for the case
   * where meshes are very large and computing indicators for @em every cell
   * is too expensive, while computing indicators for only local cells is
   * acceptable. Note that if you only ask for the indicators of a certain
   * subdomain to be computed, you must nevertheless make sure that this
   * function has access to the correct node values of @em all degrees of
   * freedom. This is since the function needs to access DoF values on
   * neighboring cells as well, even if they belong to a different subdomain.
   *
   * The @p material_id parameter has a similar meaning: if not set to its
   * default value (which is numbers::invalid_material_id), it means that
   * indicators will only be computed for cells with this particular material
   * id.
   *
   * The @p n_threads parameter used to indicate the number of threads to be
   * used to compute the error estimator. This parameter is now ignored, with
   * the number of threads determined automatically. The parameter is retained
   * for compatibility with old versions of the library.
   *
   * @note If the DoFHandler object given as an argument to this function
   * builds on a parallel::distributed::Triangulation, this function skips
   * computations on all cells that are not locally owned. In that case, the
   * only valid value for the subdomain_id argument (besides the invalid
   * value) is the subdomain id that is associated with the currently
   * processor, as reported by
   * parallel::distributed::Triangulation::locally_owned_subdomain(). Even
   * though nothing is computed on cells that we don't locally own, the error
   * indicator vector must still have a length equal to the number of active
   * cell in the mesh as reported by
   * parallel::distributed::Triangulation::n_locally_owned_active_cells().
   */
  template <typename InputVector, class DH>
  static void estimate (const Mapping<dim, spacedim>      &mapping,
                        const DH                &dof,
                        const Quadrature<dim-1> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const InputVector       &solution,
                        Vector<float>           &error,
                        const ComponentMask &component_mask = ComponentMask(),
                        const Function<spacedim>     *coefficients   = 0,
                        const unsigned int       n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id       material_id = numbers::invalid_material_id);

  /**
   * Calls the @p estimate function, see above, with
   * <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <typename InputVector, class DH>
  static void estimate (const DH                &dof,
                        const Quadrature<dim-1> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const InputVector       &solution,
                        Vector<float>           &error,
                        const ComponentMask &component_mask = ComponentMask(),
                        const Function<spacedim>     *coefficients   = 0,
                        const unsigned int       n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id       material_id = numbers::invalid_material_id);

  /**
   * Same function as above, but accepts more than one solution vector and
   * returns one error vector for each solution vector. For the reason of
   * existence of this function, see the general documentation of this class.
   *
   * Since we do not want to force the user of this function to copy around
   * their solution vectors, the vector of solution vectors takes pointers to
   * the solutions, rather than being a vector of vectors. This makes it
   * simpler to have the solution vectors somewhere in memory, rather than to
   * have them collected somewhere special. (Note that it is not possible to
   * construct of vector of references, so we had to use a vector of
   * pointers.)
   */
  template <typename InputVector, class DH>
  static void estimate (const Mapping<dim, spacedim>          &mapping,
                        const DH                    &dof,
                        const Quadrature<dim-1>     &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const std::vector<const InputVector *> &solutions,
                        std::vector<Vector<float>*> &errors,
                        const ComponentMask     &component_mask = ComponentMask(),
                        const Function<spacedim>         *coefficients   = 0,
                        const unsigned int           n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id           material_id = numbers::invalid_material_id);

  /**
   * Calls the @p estimate function, see above, with
   * <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <typename InputVector, class DH>
  static void estimate (const DH                    &dof,
                        const Quadrature<dim-1>     &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const std::vector<const InputVector *> &solutions,
                        std::vector<Vector<float>*> &errors,
                        const ComponentMask     &component_mask = ComponentMask(),
                        const Function<spacedim>         *coefficients   = 0,
                        const unsigned int           n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id          material_id = numbers::invalid_material_id);


  /**
   * Equivalent to the set of functions above, except that this one takes a
   * quadrature collection for hp finite element dof handlers.
   */
  template <typename InputVector, class DH>
  static void estimate (const Mapping<dim, spacedim>      &mapping,
                        const DH                &dof,
                        const hp::QCollection<dim-1> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const InputVector       &solution,
                        Vector<float>           &error,
                        const ComponentMask &component_mask = ComponentMask(),
                        const Function<spacedim>     *coefficients   = 0,
                        const unsigned int       n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id       material_id = numbers::invalid_material_id);


  /**
   * Equivalent to the set of functions above, except that this one takes a
   * quadrature collection for hp finite element dof handlers.
   */
  template <typename InputVector, class DH>
  static void estimate (const DH                &dof,
                        const hp::QCollection<dim-1> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const InputVector       &solution,
                        Vector<float>           &error,
                        const ComponentMask &component_mask = ComponentMask(),
                        const Function<spacedim>     *coefficients   = 0,
                        const unsigned int       n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id       material_id = numbers::invalid_material_id);


  /**
   * Equivalent to the set of functions above, except that this one takes a
   * quadrature collection for hp finite element dof handlers.
   */
  template <typename InputVector, class DH>
  static void estimate (const Mapping<dim, spacedim>          &mapping,
                        const DH                    &dof,
                        const hp::QCollection<dim-1> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const std::vector<const InputVector *> &solutions,
                        std::vector<Vector<float>*> &errors,
                        const ComponentMask     &component_mask = ComponentMask(),
                        const Function<spacedim>         *coefficients   = 0,
                        const unsigned int           n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id           material_id = numbers::invalid_material_id);


  /**
   * Equivalent to the set of functions above, except that this one takes a
   * quadrature collection for hp finite element dof handlers.
   */
  template <typename InputVector, class DH>
  static void estimate (const DH                    &dof,
                        const hp::QCollection<dim-1> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const std::vector<const InputVector *> &solutions,
                        std::vector<Vector<float>*> &errors,
                        const ComponentMask     &component_mask = ComponentMask(),
                        const Function<spacedim>    *coefficients   = 0,
                        const unsigned int           n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id           material_id = numbers::invalid_material_id);


  /**
   * Exception
   */
  DeclException0 (ExcInvalidBoundaryIndicator);
  /**
   * Exception
   */
  DeclException0 (ExcInvalidComponentMask);
  /**
   * Exception
   */
  DeclException0 (ExcInvalidCoefficient);
  /**
   * Exception
   */
  DeclException0 (ExcInvalidBoundaryFunction);
  /**
   * Exception
   */
  DeclException2 (ExcIncompatibleNumberOfElements,
                  int, int,
                  << "The number of elements " << arg1 << " and " << arg2
                  << " of the vectors do not match!");
  /**
   * Exception
   */
  DeclException0 (ExcInvalidSolutionVector);
  /**
   * Exception
   */
  DeclException0 (ExcNoSolutions);
};



/**
 * This is a specialization of the general template for 1d. The implementation
 * is sufficiently different for 1d to justify this specialization. The basic
 * difference between 1d and all other space dimensions is that in 1d, there
 * are no faces of cells, just the vertices between line segments, so we have
 * to compute the jump terms differently. However, this class offers exactly
 * the same public functions as the general template, so that a user will not
 * see any difference.
 *
 * @author Wolfgang Bangerth, 1998, 2004.
 */
template <int spacedim>
class KellyErrorEstimator<1,spacedim>
{
public:
  /**
   * Implementation of the error estimator described above. You may give a
   * coefficient, but there is a default value which denotes the constant
   * coefficient with value one. The coefficient function may either be a
   * scalar one, in which case it is used for all components of the finite
   * element, or a vector-valued one with as many components as there are in
   * the finite element; in the latter case, each component is weighted by the
   * respective component in the coefficient.
   *
   * You might give a list of components you want to evaluate, in case the
   * finite element used by the DoFHandler object is vector-valued. You then
   * have to set those entries to true in the bit-vector @p component_mask for
   * which the respective component is to be used in the error estimator. The
   * default is to use all components, which is done by either providing a
   * bit-vector with all-set entries, or an empty bit-vector. All the other
   * parameters are as in the general case used for 2d and higher.
   *
   * The parameter n_threads is no longer used and will be ignored.
   * Multithreading is not presently
   * implemented for 1d, but we retain the respective parameter for
   * compatibility with the function signature in the general case.
   */
  template <typename InputVector, class DH>
  static void estimate (const Mapping<1,spacedim>  &mapping,
                        const DH   &dof,
                        const Quadrature<0> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const InputVector       &solution,
                        Vector<float>           &error,
                        const ComponentMask &component_mask = ComponentMask(),
                        const Function<spacedim>     *coefficients   = 0,
                        const unsigned int       n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id       material_id = numbers::invalid_material_id);

  /**
   * Calls the @p estimate function, see above, with
   * <tt>mapping=MappingQ1<1>()</tt>.
   */
  template <typename InputVector, class DH>
  static void estimate (const DH   &dof,
                        const Quadrature<0> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const InputVector       &solution,
                        Vector<float>           &error,
                        const ComponentMask &component_mask = ComponentMask(),
                        const Function<spacedim>     *coefficients   = 0,
                        const unsigned int       n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id       material_id = numbers::invalid_material_id);

  /**
   * Same function as above, but accepts more than one solution vectors and
   * returns one error vector for each solution vector. For the reason of
   * existence of this function, see the general documentation of this class.
   *
   * Since we do not want to force the user of this function to copy around
   * their solution vectors, the vector of solution vectors takes pointers to
   * the solutions, rather than being a vector of vectors. This makes it
   * simpler to have the solution vectors somewhere in memory, rather than to
   * have them collected somewhere special. (Note that it is not possible to
   * construct of vector of references, so we had to use a vector of
   * pointers.)
   */
  template <typename InputVector, class DH>
  static void estimate (const Mapping<1,spacedim>          &mapping,
                        const DH       &dof,
                        const Quadrature<0>     &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const std::vector<const InputVector *> &solutions,
                        std::vector<Vector<float>*> &errors,
                        const ComponentMask     &component_mask = ComponentMask(),
                        const Function<spacedim>         *coefficients   = 0,
                        const unsigned int           n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id           material_id = numbers::invalid_material_id);

  /**
   * Calls the @p estimate function, see above, with
   * <tt>mapping=MappingQ1<1>()</tt>.
   */
  template <typename InputVector, class DH>
  static void estimate (const DH       &dof,
                        const Quadrature<0>     &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const std::vector<const InputVector *> &solutions,
                        std::vector<Vector<float>*> &errors,
                        const ComponentMask     &component_mask = ComponentMask(),
                        const Function<spacedim>         *coefficients   = 0,
                        const unsigned int           n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id           material_id = numbers::invalid_material_id);


  /**
   * Equivalent to the set of functions above, except that this one takes a
   * quadrature collection for hp finite element dof handlers.
   */
  template <typename InputVector, class DH>
  static void estimate (const Mapping<1,spacedim>      &mapping,
                        const DH                &dof,
                        const hp::QCollection<0> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const InputVector       &solution,
                        Vector<float>           &error,
                        const ComponentMask &component_mask = ComponentMask(),
                        const Function<spacedim>     *coefficients   = 0,
                        const unsigned int       n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id       material_id = numbers::invalid_material_id);


  /**
   * Equivalent to the set of functions above, except that this one takes a
   * quadrature collection for hp finite element dof handlers.
   */
  template <typename InputVector, class DH>
  static void estimate (const DH                &dof,
                        const hp::QCollection<0> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const InputVector       &solution,
                        Vector<float>           &error,
                        const ComponentMask &component_mask = ComponentMask(),
                        const Function<spacedim>     *coefficients   = 0,
                        const unsigned int       n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id       material_id = numbers::invalid_material_id);


  /**
   * Equivalent to the set of functions above, except that this one takes a
   * quadrature collection for hp finite element dof handlers.
   */
  template <typename InputVector, class DH>
  static void estimate (const Mapping<1,spacedim>          &mapping,
                        const DH                    &dof,
                        const hp::QCollection<0> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const std::vector<const InputVector *> &solutions,
                        std::vector<Vector<float>*> &errors,
                        const ComponentMask     &component_mask = ComponentMask(),
                        const Function<spacedim>         *coefficients   = 0,
                        const unsigned int           n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id           material_id = numbers::invalid_material_id);


  /**
   * Equivalent to the set of functions above, except that this one takes a
   * quadrature collection for hp finite element dof handlers.
   */
  template <typename InputVector, class DH>
  static void estimate (const DH                    &dof,
                        const hp::QCollection<0> &quadrature,
                        const typename FunctionMap<spacedim>::type &neumann_bc,
                        const std::vector<const InputVector *> &solutions,
                        std::vector<Vector<float>*> &errors,
                        const ComponentMask     &component_mask = ComponentMask(),
                        const Function<spacedim>         *coefficients   = 0,
                        const unsigned int           n_threads = numbers::invalid_unsigned_int,
                        const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
                        const types::material_id           material_id = numbers::invalid_material_id);

  /**
   * Exception
   */
  DeclException0 (ExcInvalidBoundaryIndicator);
  /**
   * Exception
   */
  DeclException0 (ExcInvalidComponentMask);
  /**
   * Exception
   */
  DeclException0 (ExcInvalidCoefficient);
  /**
   * Exception
   */
  DeclException0 (ExcInvalidBoundaryFunction);
  /**
   * Exception
   */
  DeclException2 (ExcIncompatibleNumberOfElements,
                  int, int,
                  << "The number of elements " << arg1 << " and " << arg2
                  << " of the vectors do not match!");
  /**
   * Exception
   */
  DeclException0 (ExcInvalidSolutionVector);
  /**
   * Exception
   */
  DeclException0 (ExcNoSolutions);
};



DEAL_II_NAMESPACE_CLOSE

#endif
