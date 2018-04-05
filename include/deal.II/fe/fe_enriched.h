// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#ifndef dealii_fe_enriched_h
#define dealii_fe_enriched_h

#include <deal.II/base/config.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/quadrature.h>

#include <vector>
#include <utility>
#include <numeric>

DEAL_II_NAMESPACE_OPEN

/**
 * Implementation of a partition of unity finite element method (PUM) by
 * Babuska and Melenk which enriches a standard
 * finite element with an enrichment function multiplied with another (usually
 * linear) finite element:
 * \f[
 * U(\mathbf x) = \sum_i N_i(\mathbf x) U_i + \sum_j N_j(\mathbf x) \sum_k F_k(\mathbf x) U_{jk}
 * \f]
 * where $ N_i(\mathbf x) $ and $ N_j(\mathbf x) $ are the underlying finite elements (including
 * the mapping from the isoparametric element to the real element); $ F_k(\mathbf x) $
 * are the scalar enrichment functions in real space (e.g. $ 1/r $, $ \exp(-r) $, etc);
 * $ U_i $ and $ U_{jk} $ are the standard and enriched DoFs. This allows to
 * include in the finite element space a priori knowledge about the partial
 * differential equation being solved which in turn improves the local
 * approximation properties of the spaces. This can be useful for highly oscillatory
 * solutions, problems with domain corners or on unbounded domains or sudden
 * changes of boundary conditions. PUM method uses finite element spaces which
 * satisfy the partition of unity property (e.g. FE_Q). Among other properties
 * this makes the resulting space to reproduce enrichment functions exactly.
 *
 * The simplest constructor of this class takes two finite element objects and an
 * enrichment function to be used. For example
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
 * or the finite element field requires evaluation of gradients (gradients and hessians)
 * of the enrichment functions:
 * @f{align*}{
 *   U(\mathbf x)
 *     &= \sum_i N_i(\mathbf x) U_i
 *     + \sum_{j,k} N_j(\mathbf x) F_k(\mathbf x) U_{jk} \\
 *   \mathbf \nabla U(\mathbf x)
 *     &= \sum_i \mathbf \nabla N_i(\mathbf x) U_i
 *     + \sum_{j,k} \left[\mathbf \nabla N_j(\mathbf x) F_k(\mathbf x) +
 *                        N_j(\mathbf x) \mathbf \nabla F_k(\mathbf x) \right] U_{jk} \\
 *   \mathbf \nabla \mathbf \nabla U(\mathbf x)
 *     &= \sum_i \mathbf \nabla \mathbf \nabla N_i(\mathbf x) U_i
 *     + \sum_{j,k} \left[\mathbf \nabla \mathbf \nabla N_j(\mathbf x) F_k(\mathbf x) +
 *                        \mathbf \nabla F_k(\mathbf x) \mathbf \nabla N_j(\mathbf x) +
 *                        \mathbf \nabla N_j(\mathbf x) \mathbf \nabla F_k(\mathbf x) +
 *                        N_j(\mathbf x) \mathbf \nabla \mathbf \nabla F_k(\mathbf x) \right] U_{jk}
 * @f}
 *
 * <h3>Using enriched and non-enriched FEs together</h3>
 *
 * In most applications it is beneficial to introduce enrichments only in
 * some part of the domain (e.g. around a crack tip) and use standard FE (e.g. FE_Q)
 * elsewhere.
 * This can be achieved by using the hp finite element framework in deal.II
 * that allows for the use of different elements on different cells. To make
 * the resulting space $C^0$ continuous, it is then necessary for the DoFHandler
 * class and DoFTools::make_hanging_node_constraints() function to be able to
 * figure out what to do at the interface between enriched and non-enriched
 * cells. Specifically, we want the degrees of freedom corresponding to
 * enriched shape functions to be zero at these interfaces. These classes and
 * functions can not to do this automatically, but the effect can be achieved
 * by using not just a regular FE_Q on cells without enrichment, but to wrap
 * the FE_Q into an FE_Enriched object <i>without actually enriching it</i>.
 * This can be done as follows:
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
 * @Article{Davydov2017,
 *  author    = {Denis Davydov and Tymofiy Gerasimov and Jean-Paul Pelteret and Paul Steinmann},
 *  title     = {Convergence study of the h-adaptive PUM and the hp-adaptive FEM applied to eigenvalue problems in quantum mechanics},
 *  journal   = {Advanced Modeling and Simulation in Engineering Sciences},
 *  year      = {2017},
 *  volume    = {4},
 *  number    = {1},
 *  pages     = {7},
 *  month     = {Dec},
 *  issn      = {2213-7467},
 *  day       = {12},
 *  doi       = {10.1186/s40323-017-0093-0},
 *  url       = {https://doi.org/10.1186/s40323-017-0093-0},
 * }
 * @endcode
 * The PUM was introduced in
 * @code{.bib}
 * @Article{Melenk1996,
 *   Title                    = {The partition of unity finite element method: Basic theory and applications },
 *   Author                   = {Melenk, J.M. and Babu\v{s}ka, I.},
 *   Journal                  = {Computer Methods in Applied Mechanics and Engineering},
 *   Year                     = {1996},
 *   Number                   = {1--4},
 *   Pages                    = {289 -- 314},
 *   Volume                   = {139},
 * }
 * @Article{Babuska1997,
 *   Title                    = {The partition of unity method},
 *   Author                   = {Babu\v{s}ka, I. and Melenk, J. M.},
 *   Journal                  = {International Journal for Numerical Methods in Engineering},
 *   Year                     = {1997},
 *   Number                   = {4},
 *   Pages                    = {727--758},
 *   Volume                   = {40},
 * }
 * @endcode
 *
 * <h3>Implementation</h3>
 *
 * The implementation of the class is based on FESystem which is aggregated as
 * a private member. The simplest constructor <code> FE_Enriched<dim> fe(FE_Q<dim>(2), FE_Q<dim>(1),function)</code>
 * will internally initialize FESystem as
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
template <int dim, int spacedim=dim>
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
  FE_Enriched (const FiniteElement<dim,spacedim> &fe_base,
               const FiniteElement<dim,spacedim> &fe_enriched,
               const Function<spacedim>      *enrichment_function);

  /**
   * Constructor which only wraps the base FE @p fe_base.
   * As for the enriched finite element space, FE_Nothing is used.
   * Continuity constraints will be automatically generated when
   * this non-enriched element is used in conjunction with enriched finite element
   * within the hp::DoFHandler.
   *
   * See the discussion in the class documentation on how to use this element
   * in the context of hp finite element methods.
   */
  FE_Enriched (const FiniteElement<dim,spacedim> &fe_base);

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
   * FE_Enriched<dim> fe
   * (&fe_base,
   * {&fe_1, &fe_2},
   * {{[=] (const typename Triangulation<dim>::cell_iterator &) -> const Function<dim> * {return &fe_1_function1;},
   *   [=] (const typename Triangulation<dim>::cell_iterator &) -> const Function<dim> * {return &fe_1_function2;}},
   *  {[=] (const typename Triangulation<dim>::cell_iterator &) -> const Function<dim> * {return &fe_2_function;}}});
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
  FE_Enriched (const FiniteElement<dim,spacedim> *fe_base,
               const std::vector<const FiniteElement<dim,spacedim> * > &fe_enriched,
               const std::vector<std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > > > &functions);

private:
  /**
   * The most general private constructor. The first two input parameters are
   * consistent with those in FESystem. It is used internally only with
   * <code>multiplicities[0]=1</code>, which is a logical requirement for this finite element.
   */
  FE_Enriched (const std::vector< const FiniteElement< dim, spacedim > * > &fes,
               const std::vector< unsigned int > &multiplicities,
               const std::vector<std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > > > &functions);
public:

  virtual
  std::unique_ptr<FiniteElement<dim,spacedim> >
  clone() const;

  virtual
  UpdateFlags
  requires_update_flags (const UpdateFlags update_flags) const;

  /**
   * Return a string that identifies a finite element.
   */
  virtual std::string get_name () const;

  /**
   * Access to a composing element. The index needs to be smaller than the
   * number of base elements. In the context of this class, the number of
   * base elements is always more than one: a non-enriched element plus an
   * element to be enriched, which could be FE_Nothing.
   */
  virtual const FiniteElement<dim,spacedim> &
  base_element (const unsigned int index) const;

  /**
   * Return the value of the @p ith shape function at the point @p p. @p p is a
   * point on the reference element.
   *
   * This function returns meaningful values only for non-enriched element as
   * real-space enrichment requires evaluation of the function at the point in
   * real-space.
   */
  virtual double shape_value(const unsigned int      i,
                             const Point< dim >     &p) const;

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
  get_restriction_matrix (const unsigned int child,
                          const RefinementCase<dim> &refinement_case=RefinementCase<dim>::isotropic_refinement) const;

  /**
   * Embedding matrix between grids.
   *
   * This function only makes sense when all child elements are also enriched
   * using the same function(s) as the parent element.
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix (const unsigned int child,
                           const RefinementCase<dim> &refinement_case=RefinementCase<dim>::isotropic_refinement) const;

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
  virtual bool hp_constraints_are_implemented () const;

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
  get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                 FullMatrix<double>       &matrix) const;

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
  get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;

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
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_line_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_quad_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Return whether this element dominates the one given as argument when they
   * meet at a common face, whether it is the other way around, whether
   * neither dominates, or if either could dominate.
   *
   * For a definition of domination, see FiniteElementBase::Domination and in
   * particular the
   * @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;
  //@}


  /**
   * Return enrichment functions
   */
  const std::vector<std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > > >
  get_enrichments() const;

  /**
   * Return the underlying FESystem object.
   */
  const FESystem<dim,spacedim> &
  get_fe_system() const;

protected:

  /**
   * A class to hold internal data needed for evaluation of this FE at quadrature points.
   */
  class InternalData : public FiniteElement<dim,spacedim>::InternalDataBase
  {
  public:
    /**
     * For each Finite Element (base number) and each enrichment function (base_index)
     * this struct will contain values, gradients and hessians of the enrichment functions.
     */
    struct EnrichmentValues
    {
      std::vector<double> values;
      std::vector<Tensor<1,spacedim> > gradients;
      std::vector<SymmetricTensor<2, spacedim> > hessians;
    };

    /**
     * Constructor. Is used inside setup_data to wrap FESystem's internal
     * data object. The former is called from get_data, get_subface_data and
     * get_face_data which FE_Enriched has to implement.
     *
     * Since FESystem::get_data(), FESystem::get_face_data() and FESystem::get_subface_data()
     * just create an object and return a pointer to it (i.e. they don't retain
     * ownership), we store the cast result in a std::unique_ptr to indicate
     * that InternalData owns the object.
     */
    InternalData ( std::unique_ptr<typename FESystem<dim,spacedim>::InternalData> fesystem_data);

    /**
     * Give read-access to the pointer to a @p InternalData of the @p
     * <code>base_no</code>th base element of FESystem's data.
     */
    typename FiniteElement<dim,spacedim>::InternalDataBase &
    get_fe_data (const unsigned int base_no) const;

    /**
     * Give read-access to the pointer to an object into which the
     * <code>base_no</code>th base element will write its output when calling
     * FiniteElement::fill_fe_values() and similar functions.
     */
    internal::FEValuesImplementation::FiniteElementRelatedData<dim,spacedim> &
    get_fe_output_object (const unsigned int base_no) const;

    /**
     * Aggregate FESystem's internal data. It is used every time
     * we call FESystem's fill_fe_values() and alike.
     */
    std::unique_ptr<typename FESystem<dim,spacedim>::InternalData> fesystem_data;

    /**
     * For each FE used in enrichment (base number <code>i</code>) and each enrichment function
     * (base multiplicity <code>j</code>), <code>enrichment_values[i][j]</code> will be used to store
     * possibly requested values, gradients and hessians of enrichment function <code>j</code>.
     *
     * The variable is made mutable as InternalData's provided to fill_fe_values and alike
     * are const.
     *
     * @note We do not want to store this information in the finite element object itself,
     * because this would mean that (i) only one FEValues object could use a finite element object at a time,
     * and (ii) that these objects could not be used in a multithreaded context.
     */
    mutable std::vector<std::vector<EnrichmentValues> > enrichment;
  };

  /**
   * For each finite element @p i used in enrichment and each enrichment function
   * @p j associated with it (essentially its multiplicity),
   * @p base_no_mult_local_enriched_dofs[i][j] contains the associated local DoFs
   * on the FE_Enriched finite element.
   */
  std::vector<std::vector<std::vector<unsigned int> > > base_no_mult_local_enriched_dofs;

  /**
   * Enrichment functions.
   * The size of the first vector is the same as the number of FiniteElement spaces used
   * with enrichment. Whereas the size of the inner vector corresponds to the number
   * of enrichment functions associated with a single FiniteElement.
   */
  const std::vector<std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > > > enrichments;

  /**
   * Auxiliary variable used to distinguish between the case when we do enrichment
   * and when the class simply wraps another FiniteElement.
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
  std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase>
  setup_data (std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase> &&fes_data,
              const UpdateFlags      flags,
              const Quadrature<dim_1> &quadrature) const;

  /**
   * Prepare internal data structures and fill in values independent of the
   * cell. Returns a pointer to an object of which the caller of this function
   * (FEValues) then has to assume ownership (which includes destruction when it is no
   * more needed).
   */
  virtual std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase>
  get_data (const UpdateFlags      flags,
            const Mapping<dim,spacedim>    &mapping,
            const Quadrature<dim> &quadrature,
            dealii::internal::FEValuesImplementation::FiniteElementRelatedData< dim, spacedim > &output_data) const;

  virtual std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase>
  get_face_data (const UpdateFlags      update_flags,
                 const Mapping<dim,spacedim>    &mapping,
                 const Quadrature<dim-1> &quadrature,
                 dealii::internal::FEValuesImplementation::FiniteElementRelatedData< dim, spacedim >        &output_data) const;

  virtual std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase>
  get_subface_data (const UpdateFlags      update_flags,
                    const Mapping<dim,spacedim>    &mapping,
                    const Quadrature<dim-1> &quadrature,
                    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  virtual
  void fill_fe_values (const typename Triangulation<dim, spacedim>::cell_iterator &cell,
                       const CellSimilarity::Similarity cell_similarity,
                       const Quadrature<dim> &quadrature,
                       const Mapping<dim, spacedim> &mapping,
                       const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
                       const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &mapping_data,
                       const typename FiniteElement<dim,spacedim>::InternalDataBase &fe_internal,
                       dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data
                      ) const;

  virtual
  void
  fill_fe_face_values ( const typename Triangulation<dim, spacedim>::cell_iterator   &cell,
                        const unsigned int face_no,
                        const Quadrature<dim-1> &quadrature,
                        const Mapping<dim, spacedim> &mapping,
                        const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
                        const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &mapping_data,
                        const typename FiniteElement<dim,spacedim>::InternalDataBase &fe_internal,
                        dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data
                      ) const;

  virtual
  void
  fill_fe_subface_values (const typename Triangulation< dim, spacedim >::cell_iterator &cell,
                          const unsigned int face_no,
                          const unsigned int sub_no,
                          const Quadrature<dim-1> &quadrature,
                          const Mapping<dim, spacedim> &mapping,
                          const typename Mapping< dim, spacedim >::InternalDataBase &mapping_internal,
                          const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &mapping_data,
                          const typename FiniteElement<dim,spacedim>::InternalDataBase &fe_internal,
                          dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data
                         ) const;

private:
  /**
   * This function sets up the index table for the system as well as @p
   * restriction and @p prolongation matrices.
   */
  void initialize (const std::vector<const FiniteElement<dim,spacedim>*> &fes,
                   const std::vector<unsigned int> &multiplicities);

  /**
   * The underlying FESystem object.
   */
  const std::unique_ptr<const FESystem<dim,spacedim> > fe_system;

  /**
   * After calling fill_fe_(face/subface_)values this function
   * implements the chain rule to multiply stored shape values/gradient/hessians
   * by those of enrichment function evaluated at quadrature points.
   */
  template <int dim_1>
  void
  multiply_by_enrichment (const Quadrature<dim_1>                                       &quadrature,
                          const InternalData                                            &fe_data,
                          const internal::FEValuesImplementation::MappingRelatedData<dim,spacedim>    &mapping_data,
                          const typename Triangulation< dim, spacedim >::cell_iterator  &cell,
                          internal::FEValuesImplementation::FiniteElementRelatedData<dim,spacedim>    &output_data) const;
};

//}
DEAL_II_NAMESPACE_CLOSE

#endif // dealii_fe_enriched_h
