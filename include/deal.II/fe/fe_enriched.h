// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
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

#ifndef dealii_fe_enriched_h
#define dealii_fe_enriched_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_update_flags.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>

#include <map>
#include <numeric>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Implementation of a partition of unity finite element method (PUM) by
 * Babuska and Melenk which enriches a standard
 * finite element with an enrichment function multiplied with another (usually
 * linear) finite element:
 * \f[
 * U(\mathbf x) = \sum_i N_i(\mathbf x) U_i + \sum_j N_j(\mathbf x) \sum_k
 * F_k(\mathbf x) U_{jk} \f] where $ N_i(\mathbf x) $ and $ N_j(\mathbf x) $ are
 * the underlying finite elements (including the mapping from the isoparametric
 * element to the real element); $ F_k(\mathbf x) $ are the scalar enrichment
 * functions in real space (e.g. $ 1/r $, $ \exp(-r) $, etc); $ U_i $ and $
 * U_{jk} $ are the standard and enriched DoFs. This allows to include in the
 * finite element space a priori knowledge about the partial differential
 * equation being solved which in turn improves the local approximation
 * properties of the spaces. This can be useful for highly oscillatory
 * solutions, problems with domain corners or on unbounded domains or sudden
 * changes of boundary conditions. PUM method uses finite element spaces which
 * satisfy the partition of unity property (e.g. FE_Q). Among other properties
 * this makes the resulting space to reproduce enrichment functions exactly.
 *
 * The simplest constructor of this class takes two finite element objects and
 * an enrichment function to be used. For example
 *
 * @code
 * FE_Enriched<dim> fe(FE_Q<dim>(2),
 *                     FE_Q<dim>(1),
 *                     function)
 * @endcode
 *
 * In this case, standard DoFs are distributed by <code>FE_Q<dim>(2)</code>,
 * whereas enriched DoFs are coming from a single finite element
 * <code>FE_Q<dim>(1)</code> used with a single enrichment function
 * <code>function</code>. In this case, the total number of DoFs on the
 * enriched element is the sum of DoFs from <code>FE_Q<dim>(2)</code> and
 * <code>FE_Q<dim>(1)</code>.
 *
 * As an example of an enrichment function, consider $ \exp(-x) $, which
 * leads to the following shape functions on the unit element:
 * <table>
 *   <tr>
 *     <td align="center">
 *       @image html fe_enriched_1d.png
 *     </td>
 *     <td align="center">
 *       @image html fe_enriched_h-refinement.png
 *     </td>
 *   </tr>
 *   <tr>
 *     <td align="center">
 *       1d element, base and enriched shape functions.
 *     </td>
 *     <td align="center">
 *       enriched shape function corresponding to the central vertex.
 *     </td>
 *   </tr>
 * </table>
 *
 * Note that evaluation of gradients (hessians) of the enriched shape functions
 * or the finite element field requires evaluation of gradients (gradients and
 * hessians) of the enrichment functions:
 * @f{align*}{
 *   U(\mathbf x)
 *     &= \sum_i N_i(\mathbf x) U_i
 *     + \sum_{j,k} N_j(\mathbf x) F_k(\mathbf x) U_{jk} \\
 *   \mathbf \nabla U(\mathbf x)
 *     &= \sum_i \mathbf \nabla N_i(\mathbf x) U_i
 *     + \sum_{j,k} \left[\mathbf \nabla N_j(\mathbf x) F_k(\mathbf x) +
 *                        N_j(\mathbf x) \mathbf \nabla F_k(\mathbf x) \right]
 * U_{jk} \\ \mathbf \nabla \mathbf \nabla U(\mathbf x)
 *     &= \sum_i \mathbf \nabla \mathbf \nabla N_i(\mathbf x) U_i
 *     + \sum_{j,k} \left[\mathbf \nabla \mathbf \nabla N_j(\mathbf x)
 * F_k(\mathbf x) + \mathbf \nabla F_k(\mathbf x) \mathbf \nabla N_j(\mathbf x)
 * + \mathbf \nabla N_j(\mathbf x) \mathbf \nabla F_k(\mathbf x) + N_j(\mathbf
 * x) \mathbf \nabla \mathbf \nabla F_k(\mathbf x) \right] U_{jk}
 * @f}
 *
 * <h3>Using enriched and non-enriched FEs together</h3>
 *
 * In most applications it is beneficial to introduce enrichments only in
 * some part of the domain (e.g. around a crack tip) and use standard FE (e.g.
 * FE_Q) elsewhere. This can be achieved by using the hp finite element
 * framework in deal.II that allows for the use of different elements on
 * different cells. To make the resulting space $C^0$ continuous, it is then
 * necessary for the DoFHandler class and
 * DoFTools::make_hanging_node_constraints() function to be able to figure out
 * what to do at the interface between enriched and non-enriched cells.
 * Specifically, we want the degrees of freedom corresponding to enriched shape
 * functions to be zero at these interfaces. These classes and functions can not
 * to do this automatically, but the effect can be achieved by using not just a
 * regular FE_Q on cells without enrichment, but to wrap the FE_Q into an
 * FE_Enriched object <i>without actually enriching it</i>. This can be done as
 * follows:
 * @code
 *   FE_Enriched<dim> fe_non_enriched(FE_Q<dim>(1));
 * @endcode
 * This constructor is equivalent to calling
 * @code
 *   FE_Enriched<dim> fe_non_enriched(FE_Q<dim>(1),
 *                                    FE_Nothing<dim>(1,true),
 *                                    nullptr);
 * @endcode
 * and will result in the correct constraints for enriched DoFs attributed to
 * support points on the interface between the two regions.
 *
 * <h3>References</h3>
 *
 * When using this class, please cite
 * @code{.bib}
 * @article{Davydov2017,
 *   author  = {Denis Davydov and Tymofiy Gerasimov and Jean-Paul Pelteret and
 *              Paul Steinmann},
 *   title   = {Convergence study of the h-adaptive PUM and the hp-adaptive FEM
 *              applied to eigenvalue problems in quantum mechanics},
 *   journal = {Advanced Modeling and Simulation in Engineering Sciences},
 *   year    = {2017},
 *   volume  = {4},
 *   number  = {1},
 *   pages   = {7},
 *   month   = {Dec},
 *   issn    = {2213-7467},
 *   day     = {12},
 *   doi     = {10.1186/s40323-017-0093-0},
 *   url     = {https://doi.org/10.1186/s40323-017-0093-0},
 * }
 * @endcode
 * The PUM was introduced in
 * @code{.bib}
 * @article{Melenk1996,
 *   title   = {The partition of unity finite element method: Basic theory and
 *              applications},
 *   author  = {Melenk, J.M. and Babu\v{s}ka, I.},
 *   journal = {Computer Methods in Applied Mechanics and Engineering},
 *   year    = {1996},
 *   number  = {1--4},
 *   pages   = {289 -- 314},
 *   volume  = {139},
 * }
 * @article{Babuska1997,
 *   title   = {The partition of unity method},
 *   author  = {Babu\v{s}ka, I. and Melenk, J. M.},
 *   journal = {International Journal for Numerical Methods in Engineering},
 *   year    = {1997},
 *   number  = {4},
 *   pages   = {727--758},
 *   volume  = {40},
 * }
 * @endcode
 *
 * <h3>Implementation</h3>
 *
 * The implementation of the class is based on FESystem which is aggregated as
 * a private member. The simplest constructor <code> FE_Enriched<dim>
 * fe(FE_Q<dim>(2), FE_Q<dim>(1),function)</code> will internally initialize
 * FESystem as
 *
 * @code
 * FESystem<dim> fe_system(FE_Q<dim>(2),1,
 *                         FE_Q<dim>(1),1);
 * @endcode
 *
 * Note that it would not be wise to have this class derived
 * from FESystem as the latter concatenates the given elements into different
 * components of a vector element, whereas the current class combines the given
 * elements into the same components. For instance, if two scalar elements are
 * given, the resulting element will be scalar rather than have two components
 * when doing the same with an FESystem.
 *
 * The ordering of the shape function, @p interface_constrains, the @p prolongation (embedding)
 * and the @p restriction matrices are taken from the FESystem class.
 *
 * @ingroup fe
 *
 * @author Denis Davydov, 2016.
 */
template <int dim, int spacedim = dim>
class FE_Enriched : public FiniteElement<dim, spacedim>
{
public:
  /**
   * Constructor which takes base FiniteElement @p fe_base and the enrichment
   * FiniteElement @p fe_enriched which will be multiplied by the @p enrichment_function.
   *
   * In case @p fe_enriched is other than FE_Nothing, the lifetime of the
   * @p enrichment_function must be at least as long as the FE_Enriched object.
   */
  FE_Enriched(const FiniteElement<dim, spacedim> &fe_base,
              const FiniteElement<dim, spacedim> &fe_enriched,
              const Function<spacedim> *          enrichment_function);

  /**
   * Constructor which only wraps the base FE @p fe_base.
   * As for the enriched finite element space, FE_Nothing is used.
   * Continuity constraints will be automatically generated when
   * this non-enriched element is used in conjunction with enriched finite
   * element within the hp::DoFHandler.
   *
   * See the discussion in the class documentation on how to use this element
   * in the context of hp finite element methods.
   */
  FE_Enriched(const FiniteElement<dim, spacedim> &fe_base);

  /**
   * Constructor which takes pointer to the base FiniteElement @p fe_base and
   * a vector of enriched FiniteElement's @p fe_enriched . @p fe_enriched[i]
   * finite element will be enriched with functions in @p functions[i].
   *
   * This is the most general public constructor which also allows to have
   * different enrichment functions in different disjoint parts of the domain.
   * To that end the last argument provides an association of cell iterator
   * to a Function. This is done to simplify the usage of this class when the
   * number of disjoint domains with different functions is more than a few.
   * Otherwise one would have to use different instance of this class for each
   * disjoint enriched domain.
   *
   * If you don't plan to use this feature, you can utilize C++11 lambdas to
   * define dummy functions. Below is an example which uses two functions with
   * the first element to be enriched and a single function with the second one.
   * @code
   * FE_Enriched<dim> fe(
   *   &fe_base,
   *   {&fe_1, &fe_2},
   *   {{[=] (const typename Triangulation<dim>::cell_iterator &)
   *       -> const Function<dim> *
   *     {
   *       return &fe_1_function1;
   *     },
   *     [=] (const typename Triangulation<dim>::cell_iterator &)
   *       -> const Function<dim> *
   *     {
   *       return &fe_1_function2;
   *     }},
   *    {[=] (const typename Triangulation<dim>::cell_iterator &)
   *       -> const Function<dim> *
   *     {
   *       return &fe_2_function;
   *     }
   *    }});
   * @endcode
   *
   * @note When using the same finite element for enrichment with N
   * different functions, it is advised to have the second argument of size 1
   * and the last argument of size 1 x N. The same can be achieved by providing
   * N equivalent enrichment elements while keeping the last argument of size
   * N x 1. However this will be much more computationally expensive.
   *
   * @note When using different enrichment functions on disjoint domains, no
   * checks are done by this class that the domains are actually disjoint.
   */
  FE_Enriched(
    const FiniteElement<dim, spacedim> *                     fe_base,
    const std::vector<const FiniteElement<dim, spacedim> *> &fe_enriched,
    const std::vector<std::vector<std::function<const Function<spacedim> *(
      const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
      &functions);

private:
  /**
   * The most general private constructor. The first two input parameters are
   * consistent with those in FESystem. It is used internally only with
   * <code>multiplicities[0]=1</code>, which is a logical requirement for this
   * finite element.
   */
  FE_Enriched(
    const std::vector<const FiniteElement<dim, spacedim> *> &fes,
    const std::vector<unsigned int> &                        multiplicities,
    const std::vector<std::vector<std::function<const Function<spacedim> *(
      const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
      &functions);

public:
  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * Return a string that identifies a finite element.
   */
  virtual std::string
  get_name() const override;

  /**
   * Access to a composing element. The index needs to be smaller than the
   * number of base elements. In the context of this class, the number of
   * base elements is always more than one: a non-enriched element plus an
   * element to be enriched, which could be FE_Nothing.
   */
  virtual const FiniteElement<dim, spacedim> &
  base_element(const unsigned int index) const override;

  /**
   * Return the value of the @p ith shape function at the point @p p. @p p is a
   * point on the reference element.
   *
   * This function returns meaningful values only for non-enriched element as
   * real-space enrichment requires evaluation of the function at the point in
   * real-space.
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @name Transfer matrices
   * @{
   */

  /**
   * Projection from a fine grid space onto a coarse grid space.
   *
   * This function only makes sense when all child elements are also enriched
   * using the same function(s) as the parent element.
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * Embedding matrix between grids.
   *
   * This function only makes sense when all child elements are also enriched
   * using the same function(s) as the parent element.
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  //@}

  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * Return whether this element implements hp constraints.
   *
   * This function returns @p true if and only if all its base elements return @p true
   * for this function.
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Base elements of this element will have to implement this function. They
   * may only provide interpolation matrices for certain source finite
   * elements, for example those from the same family. If they don't implement
   * interpolation from a given element, then they must throw an exception of
   * type FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented, which
   * will get propagated out from this element.
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &matrix) const override;

  /**
   * Return the matrix interpolating from a face of one element to the
   * subface of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Base elements of this element will have to implement this function. They
   * may only provide interpolation matrices for certain source finite
   * elements, for example those from the same family. If they don't implement
   * interpolation from a given element, then they must throw an exception of
   * type FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented, which
   * will get propagated out from this element.
   */
  virtual void
  get_subface_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                   const unsigned int                  subface,
                                   FullMatrix<double> &matrix) const override;

  /**
   * If, on a vertex, several finite elements are active, the hp code first
   * assigns the degrees of freedom of each of these FEs different global
   * indices. It then calls this function to find out which of them should get
   * identical values, and consequently can receive the same global DoF index.
   * This function therefore returns a list of identities between DoFs of the
   * present finite element object with the DoFs of @p fe_other, which is a
   * reference to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes which of
   * the degrees of freedom of the two finite element objects are equivalent,
   * both numbered between zero and the corresponding value of dofs_per_vertex
   * of the two finite elements. The first index of each pair denotes one of
   * the vertex dofs of the present element, whereas the second is the
   * corresponding index of the other finite element.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  //@}


  /**
   * Return enrichment functions
   */
  const std::vector<std::vector<std::function<const Function<spacedim> *(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
  get_enrichments() const;

  /**
   * Return the underlying FESystem object.
   */
  const FESystem<dim, spacedim> &
  get_fe_system() const;

protected:
  /**
   * A class to hold internal data needed for evaluation of this FE at
   * quadrature points.
   */
  class InternalData : public FiniteElement<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * For each Finite Element (base number) and each enrichment function
     * (base_index) this struct will contain values, gradients and hessians of
     * the enrichment functions.
     */
    struct EnrichmentValues
    {
      std::vector<double>                       values;
      std::vector<Tensor<1, spacedim>>          gradients;
      std::vector<SymmetricTensor<2, spacedim>> hessians;
    };

    /**
     * Constructor. Is used inside setup_data to wrap FESystem's internal
     * data object. The former is called from get_data, get_subface_data and
     * get_face_data which FE_Enriched has to implement.
     *
     * Since FESystem::get_data(), FESystem::get_face_data() and
     * FESystem::get_subface_data() just create an object and return a pointer
     * to it (i.e. they don't retain ownership), we store the cast result in a
     * std::unique_ptr to indicate that InternalData owns the object.
     */
    InternalData(std::unique_ptr<typename FESystem<dim, spacedim>::InternalData>
                   fesystem_data);

    /**
     * Give read-access to the pointer to a @p InternalData of the @p
     * <code>base_no</code>th base element of FESystem's data.
     */
    typename FiniteElement<dim, spacedim>::InternalDataBase &
    get_fe_data(const unsigned int base_no) const;

    /**
     * Give read-access to the pointer to an object into which the
     * <code>base_no</code>th base element will write its output when calling
     * FiniteElement::fill_fe_values() and similar functions.
     */
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &
    get_fe_output_object(const unsigned int base_no) const;

    /**
     * Aggregate FESystem's internal data. It is used every time
     * we call FESystem's fill_fe_values() and alike.
     */
    std::unique_ptr<typename FESystem<dim, spacedim>::InternalData>
      fesystem_data;

    /**
     * For each FE used in enrichment (base number <code>i</code>) and each
     * enrichment function (base multiplicity <code>j</code>),
     * <code>enrichment_values[i][j]</code> will be used to store possibly
     * requested values, gradients and hessians of enrichment function
     * <code>j</code>.
     *
     * The variable is made mutable as InternalData's provided to fill_fe_values
     * and alike are const.
     *
     * @note We do not want to store this information in the finite element object itself,
     * because this would mean that (i) only one FEValues object could use a
     * finite element object at a time, and (ii) that these objects could not be
     * used in a multithreaded context.
     */
    mutable std::vector<std::vector<EnrichmentValues>> enrichment;
  };

  /**
   * For each finite element @p i used in enrichment and each enrichment function
   * @p j associated with it (essentially its multiplicity),
   * @p base_no_mult_local_enriched_dofs[i][j] contains the associated local
   * DoFs on the FE_Enriched finite element.
   */
  std::vector<std::vector<std::vector<unsigned int>>>
    base_no_mult_local_enriched_dofs;

  /**
   * Enrichment functions.
   * The size of the first vector is the same as the number of FiniteElement
   * spaces used with enrichment. Whereas the size of the inner vector
   * corresponds to the number of enrichment functions associated with a single
   * FiniteElement.
   */
  const std::vector<std::vector<std::function<const Function<spacedim> *(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
    enrichments;

  /**
   * Auxiliary variable used to distinguish between the case when we do
   * enrichment and when the class simply wraps another FiniteElement.
   *
   * This variable is initialized in the constructor by looping over a vector of
   * enrichment elements and checking if all of them are FE_Nothing. If this is
   * the case, then the value is set to <code>false</code>, otherwise it is
   * <code>true</code>.
   */
  const bool is_enriched;

  /**
   * Auxiliary function called from get_data, get_face_data and
   * get_subface_data. It take internal data of FESystem object in @p fes_data
   * and the quadrature rule @p qudrature.
   *
   * This function essentially take the internal data from an instance of
   * FESystem class and wraps it into our own InternalData class which
   * additionally has objects to hold values/gradients/hessians of
   * enrichment functions at each quadrature point depending on @p flags.
   */
  template <int dim_1>
  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
  setup_data(
    std::unique_ptr<typename FESystem<dim, spacedim>::InternalData> fes_data,
    const UpdateFlags                                               flags,
    const Quadrature<dim_1> &quadrature) const;

  /**
   * Prepare internal data structures and fill in values independent of the
   * cell. Returns a pointer to an object of which the caller of this function
   * (FEValues) then has to assume ownership (which includes destruction when it
   * is no more needed).
   */
  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_face_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim - 1> &   quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_subface_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim - 1> &   quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

private:
  /**
   * This function sets up the index table for the system as well as @p
   * restriction and @p prolongation matrices.
   */
  void
  initialize(const std::vector<const FiniteElement<dim, spacedim> *> &fes,
             const std::vector<unsigned int> &multiplicities);

  /**
   * The underlying FESystem object.
   */
  const std::unique_ptr<const FESystem<dim, spacedim>> fe_system;

  /**
   * After calling fill_fe_(face/subface_)values this function
   * implements the chain rule to multiply stored shape values/gradient/hessians
   * by those of enrichment function evaluated at quadrature points.
   */
  template <int dim_1>
  void
  multiply_by_enrichment(
    const Quadrature<dim_1> &quadrature,
    const InternalData &     fe_data,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &                                                         mapping_data,
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data) const;
};



/**
 * This namespace consists of a class needed to create a collection
 * of FE_Enriched finite elements (hp::FECollection)
 * to be used with hp::DoFHandler in a domain with multiple, possibly
 * overlapping, sub-domains with individual enrichment functions.
 *
 * To create hp::FECollection, a graph coloring algorithm is used to assign
 * colors to enrichment functions before creating hp::FECollection. Hence the
 * name.
 */
namespace ColorEnriched
{
  /**
   * An alias template for predicate function which returns a
   * boolean for a Triangulation<dim,spacedim>::cell_iterator object.
   *
   * This is used by helper functions and in the implementation of
   * ColorEnriched::Helper class.
   */
  template <int dim, int spacedim = dim>
  using predicate_function = std::function<bool(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>;

#ifndef DOXYGEN
  namespace internal
  {
    /**
     * Returns true if there is a connection between subdomains in the mesh
     * associated with @p hp::DoFHandler i.e., if the subdomains share at least
     * a vertex. The two subdomains are defined by predicates provided by
     * @p predicate_1 and @p predicate_2. A predicate is a function (or
     * object of a type with an operator()) which takes in a cell iterator and
     * gives a boolean. It is said to be active in a cell if it returns true.
     *
     * An example of a custom predicate is one that checks the distance from a
     * fixed point. Note that the operator() takes in a cell iterator. Using the
     * constructor, the fixed point and the distance can be chosen.
     * @code
     * <int dim>
     * struct predicate
     * {
     *     predicate(const Point<dim> p, const int radius)
     *     :p(p),radius(radius){}
     *
     *     template <class Iterator>
     *     bool operator () (const Iterator &i)
     *     {
     *         return ( (i->center() - p).norm() < radius);
     *     }
     *
     * private:
     *     Point<dim> p;
     *     int radius;
     *
     * };
     * @endcode
     * and then the function can be used as follows to find if the subdomains
     * are connected.
     * @code
     * find_connection_between_subdomains
     * (dof_handler,
     *  predicate<dim>(Point<dim>(0,0), 1)
     *  predicate<dim>(Point<dim>(2,2), 1));
     * @endcode
     *
     * @param[in] hp::DoFHandler object
     * @param[in] predicate_1 A function (or object of a type with an
     * operator()) defining the subdomain 1. The function takes in a cell and
     * returns a boolean.
     * @param[in] predicate_2 Same as @p predicate_1 but defines subdomain 2.
     * @return A boolean "true" if the subdomains share at least a vertex.
     */
    template <int dim, int spacedim>
    bool
    find_connection_between_subdomains(
      const hp::DoFHandler<dim, spacedim> &    dof_handler,
      const predicate_function<dim, spacedim> &predicate_1,
      const predicate_function<dim, spacedim> &predicate_2);

    /**
     * Assign colors to subdomains using Graph coloring algorithm where each
     * subdomain is considered as a graph node. Subdomains which are
     * connected i.e share at least a vertex have different color. Each
     * subdomain
     * is defined using a predicate function of @p predicates.
     *
     * @param[in] dof_handler a hp::DoFHandler object
     * @param[in] predicates predicates defining the subdomains
     * @param[out] predicate_colors Colors (unsigned int) associated with each
     * subdomain.
     */
    template <int dim, int spacedim>
    unsigned int
    color_predicates(
      const hp::DoFHandler<dim, spacedim> &                 dof_handler,
      const std::vector<predicate_function<dim, spacedim>> &predicates,
      std::vector<unsigned int> &                           predicate_colors);

    /**
     * Used to construct data members @p cellwise_color_predicate_map and
     * @p fe_sets of Helper class. Inputs are hp::DoFHandler object,
     * vector of predicates and colors associated with them. Before calling
     * this function, colors can be assigned to predicates (i.e subdomains)
     * using the function color_predicates.
     *
     * Each active FE index has a set of colors associated with it.
     * A cell with an active FE index i has a set of colors given by
     * <code>fe_sets[i]</code>. An active FE index with color {a,b}
     * means that the cell has two active predicates (i.e they return true
     * for the cell) of color a and b.
     *
     * Eg: fe_sets = { {}, {1}, {2}, {1,2} } means
     * Cells with active FE index 0 have no predicates associated.
     * Cells with index 1 have a active predicate with color 1.
     * Cells with index 2 have a active predicate with color 2.
     * Cells with index 3 have active predicates with color 1 and color 2.
     *
     * A map of maps cellwise_color_predicate_map is used to associate
     * predicate colors in cells with predicate ids. For this purpose, each
     * cell is given a unique id which is stored in material id for now.
     * When the grid is refined, material id is inherited to the children, so
     * map which associates material id with color map will still be relevant.
     *
     * Now the color map can be explained with an example. If the cell with
     * material id 100 has active predicates 4 (color = 1) and 5 (color = 2),
     * the map will insert pairs (1, 4) and (2, 5) at key 100 (i.e unique id
     * of cell is mapped with a map which associates color with predicate id).
     *
     * @param[in] dof_handler hp::DoFHandler object
     * @param[in] predicates vector of predicates defining the subdomains.
     * <code>@p predicates[i]</code> returns true for a cell if it
     * belongs to subdomain with index i.
     * @param[in] predicate_colors vector of colors (unsigned int) associated
     * with each subdomain.
     * @param[out] cellwise_color_predicate_map A map of maps used to associate
     * predicate colors in cells with predicate ids.
     * @param[out] fe_sets a vector of color lists
     */
    template <int dim, int spacedim>
    void
    set_cellwise_color_set_and_fe_index(
      hp::DoFHandler<dim, spacedim> &                       dof_handler,
      const std::vector<predicate_function<dim, spacedim>> &predicates,
      const std::vector<unsigned int> &                     predicate_colors,
      std::map<unsigned int, std::map<unsigned int, unsigned int>>
        &                                  cellwise_color_predicate_map,
      std::vector<std::set<unsigned int>> &fe_sets);

    /**
     * A function that returns a vector of enrichment functions corresponding
     * to a color. The size of the vector is equal to total number of different
     * colors associated with predicates (i.e subdomains).
     *
     * Assume that a cell has a active predicates with ids 4 (color = 1) and
     * 5 (color = 2). cellwise_color_predicate_map has this information
     * provided we know the material id.
     *
     * The constructed color_enrichments is such that
     * color_enrichments[color=1](cell) will return a pointer to
     * the enrichment function with id=4, i.e. enrichments[4].
     * In other words, using the previously collected information in
     * this function we translate a vector of user provided enrichment
     * functions into a vector of functions suitable for FE_Enriched class.
     *
     * @param[in] n_colors number of colors for predicates
     * @param[in] enrichments vector of enrichment functions
     * @param[in] cellwise_color_predicate_map A map of maps used to associate
     * predicate colors in cells with predicate ids.
     * @param[out] color_enrichments A vector of functions that take in cell
     * and return a function pointer.
     */
    template <int dim, int spacedim>
    void
    make_colorwise_enrichment_functions(
      const unsigned int                                      n_colors,
      const std::vector<std::shared_ptr<Function<spacedim>>> &enrichments,
      const std::map<unsigned int, std::map<unsigned int, unsigned int>>
        &cellwise_color_predicate_map,
      std::vector<std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>>
        &color_enrichments);


    /**
     * Creates a hp::FECollection object constructed using FE_Enriched
     * elements which itself is constructed using color enrichment functions
     * and is of size equal to number of colors.
     *
     * @param[in] fe_sets a vector of color lists
     * @param[in] color_enrichments A vector of functions that take in cell
     * and return a function pointer.
     * @param[in] fe_base base FiniteElement
     * @param[in] fe_enriched enriched FiniteElements
     * @param[in] fe_nothing a finite element with zero degrees of freedom
     * @param[out] fe_collection a collection of
     * finite elements
     */
    template <int dim, int spacedim>
    void
    make_fe_collection_from_colored_enrichments(
      const unsigned int n_colors,
      const std::vector<std::set<unsigned int>>
        &fe_sets, // total list of color sets possible
      const std::vector<std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>>
        &color_enrichments, // color wise enrichment functions
      const FiniteElement<dim, spacedim> &fe_base, // basic FE element
      const FiniteElement<dim, spacedim>
        &fe_enriched, // FE multiplied by enrichment function
      const FE_Nothing<dim, spacedim> &fe_nothing,
      hp::FECollection<dim, spacedim> &fe_collection);
  }    // namespace internal
#endif // DOXYGEN



  /**
   * ColorEnriched::Helper class creates a collection of FE_Enriched finite
   * elements (hp::FECollection) to be used with hp::DoFHandler in a domain
   * with multiple, possibly overlapping, sub-domains with individual
   * enrichment functions. Note that the overlapping regions may have
   * multiple enrichment functions associated with them. This is implemented
   * using a general constructor of FE_Enriched object which allows different
   * enrichment functions.
   *
   * Consider a domain with multiple enriched sub-domains
   * which are disjoint i.e. not connected with each other.
   * To ensure $C^0$ continuity at the interface between
   * the enriched sub-domain (characterized by a single enrichment
   * function) and the non-enriched domain, we can use an FE_Enriched
   * object in the enriched sub-domain and in the non-enriched domain
   * a standard finite element (eg: FE_Q) wrapped into an FE_Enriched
   * object (which internally uses a dominating FE_Nothing object).
   * Refer to the documentation on FE_Enriched for more
   * information on this. It is to be noted that an FE_Enriched
   * object is constructed using a base FE
   * (FiniteElement objects) and one or more
   * enriched FEs. FE_Nothing is a dummy enriched FE.
   *
   * The situation becomes more
   * complicated when two enriched sub-domains
   * share an interface. When the number of enrichment functions are
   * same for the sub-domains, FE_Enriched object of one sub-domain
   * is constructed such that each enriched FE is paired (figuratively) with a
   * FE_Nothing in the FE_Enriched object of the other sub-domain.
   * For example, let the FEs fe_enr1 and fe_enr2, which will be
   * used with enrichment functions, correspond
   * to the two sub-domains. Then the FE_Enriched objects of the two
   * sub-domains are built using
   * [fe_base, fe_enr1, fe_nothing] and
   * [fe_base, fe_nothing, fe_enr2] respectively.
   * Note that the size of the vector of enriched FEs
   * (used in FE_Enriched constructor) is equal to 2, the
   * same as the number of enrichment functions. When the number of enrichment
   * functions is not the same, additional enriched FEs are paired
   * with FE_Nothing. This ensures that the enriched DOF's at the interface
   * are set to zero by the DoFTools::make_hanging_node_constraints() function.
   * Using these two strategies, we construct the appropriate FE_Enriched
   * using the general constructor. Note that this is
   * done on a mesh without hanging nodes.
   *
   * Now consider a domain with multiple sub-domains which may share
   * an interface with each other. As discussed previously,
   * the number of enriched FEs in the FE_Enriched object of each
   * sub-domain needs to be equal to the number of sub-domains. This is because
   * we are not using the information of how the domains are connected
   * and any sub-domain may share interface with any other sub-domain (not
   * considering overlaps for now!). However, in general, a given sub-domain
   * shares an interface only with a few sub-domains. This warrants
   * the use of a graph coloring algorithm to reduce
   * the size of the vector of enriched FEs
   * (used in the FE_Enriched constructor). By giving the sub-domains
   * that share no interface the same color, a single 'std::function'
   * that returns different enrichment functions for each
   * sub-domain can be constructed. Then the size of the vector of enriched
   * FEs is equal to the number of different colors
   * used for predicates (or sub-domains).
   *
   * @note The graph coloring function, SparsityTools::color_sparsity_pattern,
   * used for assigning colors to the sub-domains
   * needs MPI (use Utilities::MPI::MPI_InitFinalize to initialize MPI
   * and the necessary Zoltan setup).
   * The coloring function, based on Zoltan, is a parallel coloring
   * algorithm but is used in serial by SparsityTools::color_sparsity_pattern.
   *
   * Construction of the Helper class needs a base FiniteElement @p fe_base,
   * an enriched FiniteElement @p fe_enriched (used for all the
   * enrichment functions), a vector of predicate
   * functions (used to define sub-domains) as well as the corresponding
   * enrichment functions. The FECollection object, a collection of FE_Enriched
   * objects to be used with an hp::DoFHandler object, can be retrieved
   * using the member function build_fe_collection which also modifies the
   * active FE indices of the hp::DoFHandler object (provided as an argument
   * to the build_fe_collection function).
   *
   * <h3>Simple example</h3>
   * Consider a domain with three sub-domains defined by predicate functions.
   * Different cells are associated with FE indices as shown in the following
   * image. The three equal-sized square-shaped sub-domains 'a', 'b'
   * and 'c' can be seen. The predicates associated with these sub-domains
   * are also labeled 'a', 'b' and 'c'.
   * The subdomains 'a' and 'b' intersect with cell labeled with FE
   * index 3. The cells in 'c' are labeled with FE
   * index 1. As can be seen, connections exist between 'a' and 'b',
   * 'b' and 'c' but 'a' and 'c' are not connected.
   *
   * \htmlonly <style>div.image
   * img[src="3source_fe_indices.png"]{width:25%;}</style> \endhtmlonly
   * @image html 3source_fe_indices.png "Active FE indices"
   *
   * As discussed before, the colors of predicates are allotted using
   * the graph coloring algorithm. Each predicate is a node in the graph and if
   * two sub-domains share an interface, the corresponding predicates
   * should be given different colors.
   * (The predicate colors are different from what is shown
   * in the image. The colors in the image are as per FE indices).
   * Predicates 'a' and 'c' can be given the same color since they
   * are not connected but the color given to 'b' has to be different from
   * 'a' and 'c'.
   *
   * The name of finite element at an index (i) of @p fe_collection
   * (hp::FECollection) can be obtained by
   * <code>fe_collection[index].get_name()</code> and is
   * show in the table below. Note that all the FE_Enriched elements
   * are of the same size and FE_Nothing<2>(dominating) is used as
   * discussed before.
   *
   * <table>
   * <tr>
   * <th>Active FE index</th>
   * <th>FiniteElement name</th>
   * </tr>
   * <tr>
   * <td>0</td>
   * <td><code>FE_Enriched<2>[FE_Q<2>(2)-FE_Nothing<2>(dominating)-FE_Nothing<2>(dominating)]</code></td>
   * </tr>
   * <tr>
   * <td>1</td>
   * <td><code>FE_Enriched<2>[FE_Q<2>(2)-FE_Q<2>(1)-FE_Nothing<2>(dominating)]</code></td>
   * </tr>
   * <tr>
   * <td>2</td>
   * <td><code>FE_Enriched<2>[FE_Q<2>(2)-FE_Q<2>(1)-FE_Q<2>(1)]</code></td>
   * </tr>
   * <tr>
   * <td>3</td>
   * <td><code>FE_Enriched<2>[FE_Q<2>(2)-FE_Nothing<2>(dominating)-FE_Q<2>(1)]</code></td>
   * </tr>
   * </table>
   *
   * The internal data members used by this class need to be available when the
   * problem is solved. This can be ensured by declaring the object static,
   * which is deallocated only when the program terminates. An alternative
   * would be to use it as a data member of the containing class. Since vector
   * of predicates and enrichment functions may not be available while
   * constructing the Helper, a 'std::shared_ptr' to Helper object can be used
   * and constructed when the predicates and enrichment functions are
   * available.
   *
   * @warning The current implementation relies on assigning each cell a
   * material id, which shall not be modified after the setup
   * and h-adaptive refinement. For a given cell, the material id is used
   * to define color predicate map, which doesn't change with refinement.
   *
   * <h3>Example usage:</h3>
   * @code
   * FE_Q<dim> fe_base(2);
   * FE_Q<dim> fe_enriched(1);
   * std::vector< predicate_function<dim> > predicates;
   * std::vector< std::shared_ptr<Function<dim>> > enrichments;
   *
   * Triangulation<dim>  triangulation;
   * hp::DoFHandler<dim> dof_handler(triangulation);
   *
   * static ColorEnriched::Helper<dim> FE_helper(fe_base,
   *                                             fe_enriched,
   *                                             predicates,
   *                                             enrichments);
   * const hp::FECollection<dim>&
   * fe_collection(FE_helper.build_fe_collection(dof_handler));
   * @endcode
   *
   * @authors Nivesh Dommaraju, Denis Davydov, 2018
   */
  template <int dim, int spacedim = dim>
  struct Helper
  {
    /**
     * Constructor for Helper class.
     *
     * @param fe_base A base FiniteElement
     * @param fe_enriched An enriched FiniteElement
     * @param predicates std::vector of predicates defining the sub-domains.
     * <code>@p predicates[i]</code> returns true for a cell if it
     * belongs to a sub-domain with index (i).
     * @param enrichments std::vector of enrichment functions
     */
    Helper(const FiniteElement<dim, spacedim> &                    fe_base,
           const FiniteElement<dim, spacedim> &                    fe_enriched,
           const std::vector<predicate_function<dim, spacedim>> &  predicates,
           const std::vector<std::shared_ptr<Function<spacedim>>> &enrichments);

    /**
     * Prepares an hp::DoFHandler object. The active FE indices of
     * mesh cells are initialized to work with
     * ColorEnriched::Helper<dim,spacedim>::fe_collection.
     *
     * @param dof_handler an hp::DoFHandler object
     * @return hp::FECollection, a collection of
     * finite elements needed by @p dof_handler.
     */
    const hp::FECollection<dim, spacedim> &
    build_fe_collection(hp::DoFHandler<dim, spacedim> &dof_handler);

  private:
    /**
     * Contains a collection of FiniteElement objects needed by an
     * hp::DoFHandler object.
     */
    hp::FECollection<dim, spacedim> fe_collection;

    /**
     * A base FiniteElement used for constructing FE_Enriched
     * object required by ColorEnriched::Helper<dim,spacedim>::fe_collection.
     */
    const FiniteElement<dim, spacedim> &fe_base;

    /**
     * An enriched FiniteElement used for constructing FE_Enriched
     * object required by ColorEnriched::Helper<dim,spacedim>::fe_collection.
     */
    const FiniteElement<dim, spacedim> &fe_enriched;

    /**
     * A finite element with zero degrees of freedom used for
     * constructing FE_Enriched object required by
     * ColorEnriched::Helper<dim,spacedim>::fe_collection
     */
    const FE_Nothing<dim, spacedim> fe_nothing;

    /**
     * std::vector of predicates defining the sub-domains.
     * <code>@p predicates[i]</code> returns true for a cell if it
     * belongs to a sub-domain with index (i).
     */
    const std::vector<predicate_function<dim, spacedim>> predicates;

    /**
     * std::vector of enrichment functions corresponding
     * to the predicates. These are needed while constructing
     * ColorEnriched::Helper<dim,spacedim>::fe_collection.
     */
    const std::vector<std::shared_ptr<Function<spacedim>>> enrichments;

    /**
     * An alias template for any callable target such as functions, lambda
     * expressions, function objects that take a
     * Triangulation<dim,spacedim>::cell_iterator
     * and return a pointer to Function<dim>. This is used to define
     * Helper<dim,spacedim>::color_enrichments
     * which returns an enrichment function
     * for a cell in Triangulation<dim,spacedim>.
     */
    using cell_iterator_function = std::function<const Function<spacedim> *(
      const typename Triangulation<dim, spacedim>::cell_iterator &)>;

    /**
     * std::vector of functions that take in a cell
     * and return a function pointer. These are needed while constructing
     * fe_collection.
     *
     * color_enrichments[i](cell_iterator) returns a pointer to
     * the correct enrichment function (i.e. whose corresponding
     * predicate has the color i) for the cell.
     */
    std::vector<cell_iterator_function> color_enrichments;

    /**
     * std::vector of colors (unsigned int) associated
     * with each sub-domain. No two connected sub-domains (i.e. sub-domains that
     * share a vertex) have the same color.
     */
    std::vector<unsigned int> predicate_colors;

    /**
     * Total number of different colors in predicate_colors
     */
    unsigned int n_colors;

    /**
     * A map of maps used to associate
     * a cell with a map that in turn associates colors of active predicates in
     * the cell with corresponding predicate ids.
     */
    std::map<unsigned int, std::map<unsigned int, unsigned int>>
      cellwise_color_predicate_map;

    /**
     * A vector of different possible color sets for a given set of
     * predicates and hp::DoFHandler object
     */
    std::vector<std::set<unsigned int>> fe_sets;
  };
} // namespace ColorEnriched

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_fe_enriched_h
