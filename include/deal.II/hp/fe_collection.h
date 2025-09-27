// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_collection_h
#define dealii_fe_collection_h

#include <deal.II/base/config.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/hp/collection.h>

#include <memory>
#include <set>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
}


namespace hp
{
  /**
   * This class acts as a collection of finite element objects used in the
   * DoFHandler.
   *
   * It implements the concepts stated in the
   * @ref hpcollection
   * topic described in the doxygen documentation.
   *
   * In addition to offering access to the elements of the collection, this
   * class provides access to the maximal number of degrees of freedom per
   * vertex, line, etc, to allow allocation of as much memory as is necessary
   * in the worst case when using the finite elements associated with the
   * cells of a triangulation.
   *
   * This class has not yet been implemented for the use in the codimension
   * one case (<tt>spacedim != dim</tt>).
   *
   * @ingroup hp hpcollection
   */
  template <int dim, int spacedim = dim>
  class FECollection : public Collection<FiniteElement<dim, spacedim>>
  {
  public:
    /**
     * Whenever p-adaptivity is considered in an hp-finite element program,
     * a hierarchy of finite elements needs to be established to determine
     * succeeding finite elements for refinement and preceding ones for
     * coarsening.
     *
     * In this struct, we supply a hierarchy that is imposed on all FECollection
     * objects by default.
     */
    struct DefaultHierarchy
    {
      /**
       * Return the index succeeding @p fe_index in the @p fe_collection.
       *
       * Once the last element of the @p fe_collection is reached, there is no element on a higher level in
       * the hierarchy and thus we return the last value.
       */
      static unsigned int
      next_index(const typename hp::FECollection<dim, spacedim> &fe_collection,
                 const unsigned int                              fe_index)
      {
        return ((fe_index + 1) < fe_collection.size()) ? fe_index + 1 :
                                                         fe_index;
      }

      /**
       * Return the index preceding @p fe_index in the @p fe_collection.
       *
       * Once the first element of the @p fe_collection is reached, there is no element on a lower level in
       * the hierarchy and thus we return the first value.
       */
      static unsigned int
      previous_index(
        const typename hp::FECollection<dim, spacedim> &fe_collection,
        const unsigned int                              fe_index)
      {
        (void)fe_collection;
        return (fe_index > 0) ? fe_index - 1 : fe_index;
      }
    };

    /**
     * Default constructor. Leads to an empty collection that can later be
     * filled using push_back(). Establishes a hierarchy of finite elements
     * corresponding to their index in the collection.
     */
    FECollection();

    /**
     * Conversion constructor. This constructor creates a FECollection from a
     * single finite element. More finite element objects can be added with
     * push_back(), if desired, though it would probably be clearer to add all
     * mappings the same way.
     */
    explicit FECollection(const FiniteElement<dim, spacedim> &fe);

    /**
     * Constructor. This constructor creates a FECollection from one or
     * more finite element objects passed to the constructor. For this
     * call to be valid, all arguments need to be of types derived
     * from class FiniteElement<dim,spacedim>.
     */
    template <class... FETypes>
    explicit FECollection(const FETypes &...fes);

    /**
     * Constructor. Same as above but for any number of elements. Pointers to
     * the elements are passed in a vector to this constructor. As above, the
     * finite element objects pointed to by the argument are not actually used
     * other than to create copies internally. Consequently, you can delete
     * these pointers immediately again after calling this constructor.
     */
    FECollection(const std::vector<const FiniteElement<dim, spacedim> *> &fes);

    /**
     * Copy constructor.
     */
    FECollection(const FECollection<dim, spacedim> &) = default;

    /**
     * Move constructor.
     *
     * @note The implementation of standard datatypes may change with different
     * libraries, so their move members may or may not be flagged non-throwing.
     * We need to explicitly set the noexcept specifier according to its
     * member variables to still get the performance benefits (and to satisfy
     * clang-tidy).
     */
    FECollection(FECollection<dim, spacedim> &&) noexcept(
      std::is_nothrow_move_constructible_v<
        std::vector<std::shared_ptr<const FiniteElement<dim, spacedim>>>>
        &&std::is_nothrow_move_constructible_v<std::function<
          unsigned int(const typename hp::FECollection<dim, spacedim> &,
                       const unsigned int)>>) = default;

    /**
     * Move assignment operator.
     */
    FECollection<dim, spacedim> &
    operator=(FECollection<dim, spacedim> &&) = default; // NOLINT

    /**
     * Equality comparison operator. All stored FiniteElement objects are
     * compared in order.
     */
    bool
    operator==(const FECollection<dim, spacedim> &fe_collection) const;

    /**
     * Non-equality comparison operator. All stored FiniteElement objects are
     * compared in order.
     */
    bool
    operator!=(const FECollection<dim, spacedim> &fe_collection) const;

    /**
     * Add a finite element. This function generates a copy of the given
     * element, i.e. you can do things like <tt>push_back(FE_Q<dim>(1));</tt>.
     * The internal copy is later destroyed by this object upon destruction of
     * the entire collection.
     *
     * When a new element is added, it needs to have the same number of vector
     * components as all other elements already in the collection.
     */
    void
    push_back(const FiniteElement<dim, spacedim> &new_fe);

    /**
     * @name Querying information about the elements in the collection
     */

    /**
     * @{
     */

    /**
     * Return the number of vector components of the finite elements in this
     * collection.  This number must be the same for all elements in the
     * collection.
     *
     * This function calls FiniteElement::n_components.  See
     * @ref GlossComponent "the glossary"
     * for more information.
     */
    unsigned int
    n_components() const;

    /**
     * Return the number of vector blocks of the finite elements in this
     * collection. While this class ensures that all elements stored in it
     * have the same number of vector components, there is no such guarantees
     * for the number of blocks each element is made up of (an element may
     * have fewer blocks than vector components; see
     * @ref GlossBlock "the glossary"
     * for more information). For example, you may have an FECollection object
     * that stores one copy of an FESystem with <code>dim</code> FE_Q objects
     * and one copy of an FE_RaviartThomas element. Both have <code>dim</code>
     * vector components but while the former has <code>dim</code> blocks the
     * latter has only one. Consequently, this function will throw an
     * assertion if the number of blocks is not the same for all elements. If
     * they are the same, this function returns the result of
     * FiniteElement::n_blocks().
     */
    unsigned int
    n_blocks() const;

    /**
     * Return the maximum of values returned by FiniteElement::get_degree()
     * over all elements of this collection.
     */
    unsigned int
    max_degree() const;

    /**
     * Return the maximal number of degrees of freedom per vertex over all
     * elements of this collection.
     */
    unsigned int
    max_dofs_per_vertex() const;

    /**
     * Return the maximal number of degrees of freedom per line over all
     * elements of this collection.
     */
    unsigned int
    max_dofs_per_line() const;

    /**
     * Return the maximal number of degrees of freedom per quad over all
     * elements of this collection.
     */
    unsigned int
    max_dofs_per_quad() const;

    /**
     * Return the maximal number of degrees of freedom per hex over all
     * elements of this collection.
     */
    unsigned int
    max_dofs_per_hex() const;

    /**
     * Return the maximal number of degrees of freedom per face over all
     * elements of this collection.
     */
    unsigned int
    max_dofs_per_face() const;

    /**
     * Return the maximal number of degrees of freedom per cell over all
     * elements of this collection.
     */
    unsigned int
    max_dofs_per_cell() const;

    /**
     * Return a mapping collection that consists of the default linear mappings
     * matching the reference cells for each hp index. More details may be found
     * in the documentation for ReferenceCell::get_default_linear_mapping().
     *
     * @note This FECollection object must remain in scope for as long as the
     * reference cell default linear mapping is in use.
     */
    const MappingCollection<dim, spacedim> &
    get_reference_cell_default_linear_mapping() const;

    /**
     * @}
     */

    /**
     * @name Functions to support hp-adaptivity
     */

    /**
     * @{
     */

    /**
     * Return whether all elements in this collection implement the hanging
     * node constraints in the new way, which has to be used to make elements
     * "hp-compatible". If this is not the case, the function returns false,
     * which implies, that at least one element in the FECollection does not
     * support the new face interface constraints. On the other hand, if this
     * method does return true, this does not imply that the hp-method will
     * work!
     *
     * This behavior is related to the fact, that FiniteElement classes,
     * which provide the new style hanging node constraints might still not
     * provide them for all possible cases. If FE_Q and FE_RaviartThomas
     * elements are included in the FECollection and both properly implement
     * the get_face_interpolation_matrix method, this method will return true.
     * But the get_face_interpolation_matrix might still fail to find an
     * interpolation matrix between these two elements.
     */
    bool
    hp_constraints_are_implemented() const;

    /**
     * This function combines the functionality of the
     * FiniteElement::hp_vertex_dof_identities() into multi-way comparisons.
     * Given a set of elements (whose indices are provided as argument), this
     * function determines identities between degrees of freedom of these
     * elements at a vertex.
     *
     * The function returns a vector of such identities, where each element of
     * the vector is a set of pairs `(fe_index,dof_index)` that identifies
     * the `fe_index` (an element of the `fes` argument to this function) of
     * an element and the `dof_index` indicates the how-manyth degree of freedom
     * of that element on a vertex participates in this identity. Now,
     * every `fe_index` can appear only once in these sets (for each identity,
     * only one degree of freedom of a finite element can be involved --
     * otherwise we would have identities between different DoFs of the same
     * element, which would make the element not unisolvent), and as a
     * consequence the function does not actually return a set of
     * `(fe_index,dof_index)` pairs for each identity, but instead a `std::map`
     * from `fe_index` to `dof_index`, which is conceptually of course
     * equivalent to a `std::set` of pairs, but in practice is easier to query.
     */
    std::vector<std::map<unsigned int, unsigned int>>
    hp_vertex_dof_identities(const std::set<unsigned int> &fes) const;

    /**
     * Same as hp_vertex_dof_indices(), except that the function treats degrees
     * of freedom on lines.
     */
    std::vector<std::map<unsigned int, unsigned int>>
    hp_line_dof_identities(const std::set<unsigned int> &fes) const;

    /**
     * Same as hp_vertex_dof_indices(), except that the function treats degrees
     * of freedom on quads.
     */
    std::vector<std::map<unsigned int, unsigned int>>
    hp_quad_dof_identities(const std::set<unsigned int> &fes,
                           const unsigned int            face_no = 0) const;


    /**
     * Return the indices of finite elements in this FECollection that dominate
     * all elements associated with the provided set of indices @p fes.
     *
     * You may find information about the domination behavior of finite elements
     * in their respective class documentation or in the implementation of their
     * inherited member function FiniteElement::compare_for_domination().
     * Consider that a finite element may or may not dominate itself (e.g.
     * FE_Nothing elements).
     *
     * For example, if a FECollection consists of
     * `{FE_Q(1),FE_Q(2),FE_Q(3),FE_Q(4)}` elements and we are looking for the
     * finite elements that dominate the middle elements of this
     * collection (i.e., @p fes is `{1,2}`), then the answer is `{FE_Q(1),FE_Q(2)`
     * and therefore this function will return their indices in the
     * FECollection, namely `{0,1}`.
     *
     * The @p codim parameter describes the codimension of the investigated
     * subspace and specifies that it is subject to this comparison. See
     * FiniteElement::compare_for_domination() for more information.
     */
    std::set<unsigned int>
    find_common_fes(const std::set<unsigned int> &fes,
                    const unsigned int            codim = 0) const;

    /**
     * Return the indices of finite elements in this FECollection that are
     * dominated by all elements associated with the provided set of indices @p fes.
     *
     * You may find information about the domination behavior of finite elements
     * in their respective class documentation or in the implementation of their
     * inherited member function FiniteElement::compare_for_domination().
     * Consider that a finite element may or may not dominate itself (e.g.
     * FE_Nothing elements).
     *
     * For example, if a FECollection consists of
     * `{FE_Q(1),FE_Q(2),FE_Q(3),FE_Q(4)}` elements and we are looking for the
     * finite elements that are dominated by the middle elements of this
     * collection (i.e., @p fes is `{1,2}`), then the answer is `{FE_Q(3),FE_Q(4)`
     * and therefore this function will return their indices in the
     * FECollection, namely `{2,3}`.
     *
     * The @p codim parameter describes the codimension of the investigated
     * subspace and specifies that it is subject to this comparison. See
     * FiniteElement::compare_for_domination() for more information.
     */
    std::set<unsigned int>
    find_enclosing_fes(const std::set<unsigned int> &fes,
                       const unsigned int            codim = 0) const;

    /**
     * Return the index of a finite element from the provided set of indices @p fes
     * that dominates all other elements of this very set.
     *
     * You may find information about the domination behavior of finite elements
     * in their respective class documentation or in the implementation of their
     * inherited member function FiniteElement::compare_for_domination().
     * Consider that a finite element may or may not dominate itself (e.g.
     * FE_Nothing elements).
     *
     * If this set consists of exactly one element, we consider it to be
     * the dominating one and return its corresponding index. Further, if the
     * function is not able to find a finite element at all, it returns
     * numbers::invalid_fe_index.
     *
     * For example, if a FECollection consists of
     * `{FE_Q(1),FE_Q(2),FE_Q(3),FE_Q(4)}` elements and we are looking for the
     * dominating finite element among the middle elements of this
     * collection (i.e., @p fes is `{1,2}`), then the answer is FE_Q(2)
     * and therefore this function will return its index in the
     * FECollection, namely `1`.
     *
     * It is of course possible that there is more than one element that
     * dominates all selected elements. For example, if the collection consists
     * of `{FE_Q(1),FE_Q(1),FE_Q(2),FE_Q(2)}` and `fes` covers all indices,
     * then one could return zero or one.  In that case, the function returns
     * either `0` or `1` since there is no tie-breaker between the two.
     *
     * The @p codim parameter describes the codimension of the investigated
     * subspace and specifies that it is subject to this comparison. See
     * FiniteElement::compare_for_domination() for more information.
     */
    unsigned int
    find_dominating_fe(const std::set<unsigned int> &fes,
                       const unsigned int            codim = 0) const;

    /**
     * Return the index of a finite element from the provided set of indices @p fes
     * that is dominated by all other elements of this very set.
     *
     * You may find information about the domination behavior of finite elements
     * in their respective class documentation or in the implementation of their
     * inherited member function FiniteElement::compare_for_domination().
     * Consider that a finite element may or may not dominate itself (e.g.
     * FE_Nothing elements).
     *
     * If this set consists of exactly one element, we consider it to be
     * the dominated one and return its corresponding index. Further, if the
     * function is not able to find a finite element at all, it returns
     * numbers::invalid_fe_index.
     *
     * For example, if a FECollection consists of
     * `{FE_Q(1),FE_Q(2),FE_Q(3),FE_Q(4)}` elements and we are looking for the
     * dominated finite element among the middle elements of this
     * collection (i.e., @p fes is `{1,2}`), then the answer is FE_Q(3)
     * and therefore this function will return its index in the
     * FECollection, namely `2`.
     *
     * It is of course possible that there is more than one element that is
     * dominated by all selected elements. For example, if the collection
     * consists of `{FE_Q(1),FE_Q(1),FE_Q(2),FE_Q(2)}` and `fes` covers all
     * indices, then one could return two or three.  In that case, the function
     * returns either `2` or `3` since there is no tie-breaker between the two.
     *
     * The @p codim parameter describes the codimension of the investigated
     * subspace and specifies that it is subject to this comparison. See
     * FiniteElement::compare_for_domination() for more information.
     */
    unsigned int
    find_dominated_fe(const std::set<unsigned int> &fes,
                      const unsigned int            codim = 0) const;

    /**
     * Return the index of a finite element from the provided set of indices @p fes
     * that dominates all other elements of this very set. If we do not succeed,
     * we extend our search on the whole collection by picking the least
     * dominating one, which is the element that describes the largest finite
     * element space of which all of the finite elements of the
     * provided set @p fes are part of.
     *
     * You may find information about the domination behavior of finite elements
     * in their respective class documentation or in the implementation of their
     * inherited member function FiniteElement::compare_for_domination().
     * Consider that a finite element may or may not dominate itself (e.g.
     * FE_Nothing elements).
     *
     * If this set consists of exactly one element, we consider it to be
     * the dominated one and return its corresponding index. Further, if the
     * function is not able to find a finite element at all, it returns
     * numbers::invalid_fe_index.
     *
     * The @p codim parameter describes the codimension of the investigated
     * subspace and specifies that it is subject to this comparison. See
     * FiniteElement::compare_for_domination() for more information.
     */
    unsigned int
    find_dominating_fe_extended(const std::set<unsigned int> &fes,
                                const unsigned int            codim = 0) const;

    /**
     * Return the index of a finite element from the provided set of indices @p fes
     * that is dominated by all other elements of this very set. If we do not
     * succeed, we extend our search on the whole collection by picking the most
     * dominated one, which is the element that describes the smallest finite
     * element space which includes all finite elements of the provided set @p fes.
     *
     * You may find information about the domination behavior of finite elements
     * in their respective class documentation or in the implementation of their
     * inherited member function FiniteElement::compare_for_domination().
     * Consider that a finite element may or may not dominate itself (e.g.
     * FE_Nothing elements).
     *
     * If this set consists of exactly one element, we consider it to be
     * the dominating one and return its corresponding index. Further, if the
     * function is not able to find a finite element at all, it returns
     * numbers::invalid_fe_index.
     *
     * The @p codim parameter describes the codimension of the investigated
     * subspace and specifies that it is subject to this comparison. See
     * FiniteElement::compare_for_domination() for more information.
     */
    unsigned int
    find_dominated_fe_extended(const std::set<unsigned int> &fes,
                               const unsigned int            codim = 0) const;

    /**
     * @}
     */

    /**
     * @name Describing hierarchical relationships between elements
     */

    /**
     * @{
     */

    /**
     * Set functions determining the hierarchy of finite elements, i.e. a
     * function @p next that returns the index of the finite element following
     * the given one, and a function @p prev returning the preceding one.
     *
     * Both functions expect an hp::FECollection to be passed along with a
     * finite element index, on whose basis the new index will be found and
     * returned.
     *
     * @note Both passed and returned indices have to be valid within the index
     * range of this collection, i.e. within [0, size()).
     */
    void
    set_hierarchy(const std::function<unsigned int(
                    const typename hp::FECollection<dim, spacedim> &,
                    const unsigned int)> &next,
                  const std::function<unsigned int(
                    const typename hp::FECollection<dim, spacedim> &,
                    const unsigned int)> &prev);

    /**
     * Set the default hierarchy corresponding to the index of each finite
     * element in the collection.
     *
     * This default hierarchy is established with functions
     * DefaultHierarchy::next_index() and DefaultHierarchy::previous_index().
     */
    void
    set_default_hierarchy();

    /**
     * Returns a sequence of FE indices that corresponds to the registered
     * hierarchy in ascending order, i.e., FE indices are sorted from lowest to
     * highest level.
     *
     * Multiple sequences of FE indices are possible with a single custom
     * hierarchy that can be registered with set_hierarchy(). This function
     * will return the sequence that contains the user-provided index
     * @p fe_index which could be located anywhere inside the sequence. The
     * default hierarchy set via set_default_hierarchy(), which corresponds to
     * FE indices in ascending order, consists of only one sequence.
     *
     * This function can be used, for example, to verify that your provided
     * hierarchy covers all elements in the desired order.
     *
     * Only one sequence of FE indices exists if the size of the returned
     * container equals the number of elements of this object, i.e.,
     * FECollection::size().
     */
    std::vector<unsigned int>
    get_hierarchy_sequence(const unsigned int fe_index = 0) const;

    /**
     * %Function returning the index of the finite element following the given
     * @p fe_index in hierarchy.
     *
     * By default, the index succeeding @p fe_index will be returned. If @p fe_index
     * already corresponds to the last index, the last index will be returned.
     * A custom hierarchy can be supplied via the member function
     * set_hierarchy().
     */
    unsigned int
    next_in_hierarchy(const unsigned int fe_index) const;

    /**
     * %Function returning the index of the finite element preceding the given
     * @p fe_index in hierarchy.
     *
     * By default, the index preceding @p fe_index will be returned. If @p fe_index
     * already corresponds to the first index, the first index will be returned.
     * A custom hierarchy can be supplied via the member function
     * set_hierarchy().
     */
    unsigned int
    previous_in_hierarchy(const unsigned int fe_index) const;

    /**
     * @}
     */

    /**
     * @name Components and blocks of elements
     */

    /**
     * @{
     */

    /**
     * Return a component mask with as many elements as this object has vector
     * components and of which exactly the one component is true that
     * corresponds to the given argument.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments. It verifies
     * that it gets the same result from every one of the elements that are
     * stored in this FECollection. If this is not the case, it throws an
     * exception.
     *
     * @param scalar An object that represents a single scalar vector
     * component of this finite element.
     * @return A component mask that is false in all components except for the
     * one that corresponds to the argument.
     */
    ComponentMask
    component_mask(const FEValuesExtractors::Scalar &scalar) const;

    /**
     * Return a component mask with as many elements as this object has vector
     * components and of which exactly the <code>dim</code> components are
     * true that correspond to the given argument.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments. It verifies
     * that it gets the same result from every one of the elements that are
     * stored in this FECollection. If this is not the case, it throws an
     * exception.
     *
     * @param vector An object that represents dim vector components of this
     * finite element.
     * @return A component mask that is false in all components except for the
     * ones that corresponds to the argument.
     */
    ComponentMask
    component_mask(const FEValuesExtractors::Vector &vector) const;

    /**
     * Return a component mask with as many elements as this object has vector
     * components and of which exactly the <code>dim*(dim+1)/2</code>
     * components are true that correspond to the given argument.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments. It verifies
     * that it gets the same result from every one of the elements that are
     * stored in this FECollection. If this is not the case, it throws an
     * exception.
     *
     * @param sym_tensor An object that represents dim*(dim+1)/2 components of
     * this finite element that are jointly to be interpreted as forming a
     * symmetric tensor.
     * @return A component mask that is false in all components except for the
     * ones that corresponds to the argument.
     */
    ComponentMask
    component_mask(
      const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

    /**
     * Given a block mask (see
     * @ref GlossBlockMask "this glossary entry"
     * ), produce a component mask (see
     * @ref GlossComponentMask "this glossary entry"
     * ) that represents the components that correspond to the blocks selected
     * in the input argument. This is essentially a conversion operator from
     * BlockMask to ComponentMask.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments. It verifies
     * that it gets the same result from every one of the elements that are
     * stored in this FECollection. If this is not the case, it throws an
     * exception.
     *
     * @param block_mask The mask that selects individual blocks of the finite
     * element
     * @return A mask that selects those components corresponding to the
     * selected blocks of the input argument.
     */
    ComponentMask
    component_mask(const BlockMask &block_mask) const;

    /**
     * Return a block mask with as many elements as this object has blocks and
     * of which exactly the one component is true that corresponds to the
     * given argument. See
     * @ref GlossBlockMask "the glossary"
     * for more information.
     *
     * @note This function will only succeed if the scalar referenced by the
     * argument encompasses a complete block. In other words, if, for example,
     * you pass an extractor for the single $x$ velocity and this object
     * represents an FE_RaviartThomas object, then the single scalar object
     * you selected is part of a larger block and consequently there is no
     * block mask that would represent it. The function will then produce an
     * exception.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments. It verifies
     * that it gets the same result from every one of the elements that are
     * stored in this FECollection. If this is not the case, it throws an
     * exception.
     *
     * @param scalar An object that represents a single scalar vector
     * component of this finite element.
     * @return A component mask that is false in all components except for the
     * one that corresponds to the argument.
     */
    BlockMask
    block_mask(const FEValuesExtractors::Scalar &scalar) const;

    /**
     * Return a component mask with as many elements as this object has vector
     * components and of which exactly the <code>dim</code> components are
     * true that correspond to the given argument. See
     * @ref GlossBlockMask "the glossary"
     * for more information.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments. It verifies
     * that it gets the same result from every one of the elements that are
     * stored in this FECollection. If this is not the case, it throws an
     * exception.
     *
     * @note The same caveat applies as to the version of the function above:
     * The extractor object passed as argument must be so that it corresponds
     * to full blocks and does not split blocks of this element.
     *
     * @param vector An object that represents dim vector components of this
     * finite element.
     * @return A component mask that is false in all components except for the
     * ones that corresponds to the argument.
     */
    BlockMask
    block_mask(const FEValuesExtractors::Vector &vector) const;

    /**
     * Return a component mask with as many elements as this object has vector
     * components and of which exactly the <code>dim*(dim+1)/2</code>
     * components are true that correspond to the given argument. See
     * @ref GlossBlockMask "the glossary"
     * for more information.
     *
     * @note The same caveat applies as to the version of the function above:
     * The extractor object passed as argument must be so that it corresponds
     * to full blocks and does not split blocks of this element.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments. It verifies
     * that it gets the same result from every one of the elements that are
     * stored in this FECollection. If this is not the case, it throws an
     * exception.
     *
     * @param sym_tensor An object that represents dim*(dim+1)/2 components of
     * this finite element that are jointly to be interpreted as forming a
     * symmetric tensor.
     * @return A component mask that is false in all components except for the
     * ones that corresponds to the argument.
     */
    BlockMask
    block_mask(const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

    /**
     * Given a component mask (see
     * @ref GlossComponentMask "this glossary entry"
     * ), produce a block mask (see
     * @ref GlossBlockMask "this glossary entry"
     * ) that represents the blocks that correspond to the components selected
     * in the input argument. This is essentially a conversion operator from
     * ComponentMask to BlockMask.
     *
     * @note This function will only succeed if the components referenced by
     * the argument encompasses complete blocks. In other words, if, for
     * example, you pass an component mask for the single $x$ velocity and
     * this object represents an FE_RaviartThomas object, then the single
     * component you selected is part of a larger block and consequently there
     * is no block mask that would represent it. The function will then
     * produce an exception.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments. It verifies
     * that it gets the same result from every one of the elements that are
     * stored in this FECollection. If this is not the case, it throws an
     * exception.
     *
     * @param component_mask The mask that selects individual components of
     * the finite element
     * @return A mask that selects those blocks corresponding to the selected
     * blocks of the input argument.
     */
    BlockMask
    block_mask(const ComponentMask &component_mask) const;

    /**
     * @}
     */

    /**
     * @name Exceptions
     * @{
     */

    /**
     * Exception
     *
     * @ingroup Exceptions
     */
    DeclException0(ExcNoFiniteElements);

    /**
     * @}
     */

  private:
    /**
     * A linear mapping collection for all reference cell types of each index
     * of this object.
     */
    std::shared_ptr<MappingCollection<dim, spacedim>>
      reference_cell_default_linear_mapping;

    /**
     * %Function returning the index of the finite element following the given
     * one in hierarchy.
     */
    std::function<unsigned int(const typename hp::FECollection<dim, spacedim> &,
                               const unsigned int)>
      hierarchy_next;

    /**
     * %Function returning the index of the finite element preceding the given
     * one in hierarchy.
     */
    std::function<unsigned int(const typename hp::FECollection<dim, spacedim> &,
                               const unsigned int)>
      hierarchy_prev;
  };



  /* --------------- inline functions ------------------- */

  template <int dim, int spacedim>
  template <class... FETypes>
  FECollection<dim, spacedim>::FECollection(const FETypes &...fes)
  {
    static_assert(
      is_base_of_all<FiniteElement<dim, spacedim>, FETypes...>::value,
      "Not all of the input arguments of this function "
      "are derived from FiniteElement<dim, spacedim>!");

    // loop over all of the given arguments and add the finite elements to
    // this collection. Inlining the definition of fe_pointers causes internal
    // compiler errors on GCC 7.1.1 so we define it separately:
    const auto fe_pointers = {
      (static_cast<const FiniteElement<dim, spacedim> *>(&fes))...};
    for (const auto p : fe_pointers)
      push_back(*p);
  }



  template <int dim, int spacedim>
  inline unsigned int
  FECollection<dim, spacedim>::n_components() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    // note that there is no need
    // here to enforce that indeed
    // all elements have the same
    // number of components since we
    // have already done this when
    // adding a new element to the
    // collection.

    return this->operator[](0).n_components();
  }



  template <int dim, int spacedim>
  inline bool
  FECollection<dim, spacedim>::operator==(
    const FECollection<dim, spacedim> &fe_collection) const
  {
    const unsigned int n_elements = this->size();
    if (n_elements != fe_collection.size())
      return false;

    for (unsigned int i = 0; i < n_elements; ++i)
      if (!(this->operator[](i) == fe_collection[i]))
        return false;

    return true;
  }



  template <int dim, int spacedim>
  inline bool
  FECollection<dim, spacedim>::operator!=(
    const FECollection<dim, spacedim> &fe_collection) const
  {
    return !(*this == fe_collection);
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_degree() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).degree);

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_vertex() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).n_dofs_per_vertex());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_line() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).n_dofs_per_line());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_quad() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).max_dofs_per_quad());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_hex() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).n_dofs_per_hex());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_face() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).max_dofs_per_face());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_cell() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).n_dofs_per_cell());

    return max;
  }


  template <int dim, int spacedim>
  bool
  FECollection<dim, spacedim>::hp_constraints_are_implemented() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    for (unsigned int i = 0; i < this->size(); ++i)
      if (this->operator[](i).hp_constraints_are_implemented() == false)
        return false;

    return true;
  }


} // namespace hp

DEAL_II_NAMESPACE_CLOSE

#endif
