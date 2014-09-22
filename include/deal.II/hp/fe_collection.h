// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

#ifndef __deal2__fe_collection_h
#define __deal2__fe_collection_h

#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/component_mask.h>

DEAL_II_NAMESPACE_OPEN

namespace hp
{

  /**
   * This class acts as a collection of finite element objects used in the
   * hp::DoFHandler. It is thus to a hp::DoFHandler what a
   * FiniteElement is to a ::DoFHandler.
   *
   * It implements the concepts stated in the @ref hpcollection module described
   * in the doxygen documentation.
   *
   * In addition to offering access to the elements of the collection, this
   * class provides access to the maximal number of degrees of freedom per
   * vertex, line, etc, to allow allocation of as much memory as is necessary in
   * the worst case when using the finite elements associated with the cells of
   * a triangulation.
   *
   * This class has not yet been implemented for the use in the codimension
   * one case (<tt>spacedim != dim </tt>).
   *
   * @ingroup hp hpcollection
   *
   * @author Wolfgang Bangerth, 2003
   */
  template <int dim, int spacedim=dim>
  class FECollection : public Subscriptor
  {
  public:
    /**
     * Default constructor. Leads
     * to an empty collection that
     * can later be filled using
     * push_back().
     */
    FECollection ();

    /**
     * Conversion constructor. This
     * constructor creates a
     * FECollection from a single
     * finite element. More finite
     * element objects can be added
     * with push_back(), if
     * desired, though it would
     * probably be clearer to add
     * all mappings the same way.
     */
    explicit FECollection (const FiniteElement<dim,spacedim> &fe);

    /**
     * Copy constructor.
     */
    FECollection (const FECollection<dim,spacedim> &fe_collection);

    /**
     * Add a finite element. This
     * function generates a copy of
     * the given element, i.e. you
     * can do things like
     * <tt>push_back(FE_Q<dim>(1));</tt>. The
     * internal copy is later
     * destroyed by this object
     * upon destruction of the
     * entire collection.
     *
     * When a new element is added,
     * it needs to have the same
     * number of vector components
     * as all other elements
     * already in the collection.
     */
    void push_back (const FiniteElement<dim,spacedim> &new_fe);

    /**
     * Get a reference to the given
     * element in this collection.
     *
     * @pre @p index must be
     * between zero and the number
     * of elements of the
     * collection.
     */
    const FiniteElement<dim,spacedim> &
    operator[] (const unsigned int index) const;

    /**
     * Return the number of finite
     * element objects stored in
     * this collection.
     */
    unsigned int size () const;

    /**
     * Return the number of vector
     * components of the finite elements in
     * this collection.  This number must
     * be the same for all elements in the
     * collection.
    *
    * This function calls
    * FiniteElement::n_components.  See
    * @ref GlossComponent "the glossary"
    * for more information.
     */
    unsigned int n_components () const;

    /**
     * Return the number of vector
     * blocks of the finite elements in
     * this collection. While this class ensures that all
     * elements stored in it have the same number of vector
     * components, there is no such guarantees for the number
     * of blocks each element is made up of (an element may
     * have fewer blocks than vector components; see
     * @ref GlossBlock "the glossary" for more information).
     * For example, you may have an FECollection object that stores
     * one copy of an FESystem with <code>dim</code> FE_Q objects
     * and one copy of an FE_RaviartThomas element. Both have
     * <code>dim</code> vector components but while the former has
     * <code>dim</code> blocks the latter has only one.
     * Consequently, this function will throw an assertion
     * if the number of blocks is not the same for all elements.
     * If they are the same, this function returns the result
     * of FiniteElement::n_blocks().
     */
    unsigned int n_blocks () const;

    /**
     * Return the maximal number of degrees
     * of freedom per vertex over all
     * elements of this collection.
     */
    unsigned int max_dofs_per_vertex () const;

    /**
     * Return the maximal number of degrees
     * of freedom per line over all elements
     * of this collection.
     */
    unsigned int max_dofs_per_line () const;

    /**
     * Return the maximal number of degrees
     * of freedom per quad over all elements
     * of this collection.
     */
    unsigned int max_dofs_per_quad () const;

    /**
     * Return the maximal number of degrees
     * of freedom per hex over all elements
     * of this collection.
     */
    unsigned int max_dofs_per_hex () const;

    /**
     * Return the maximal number of degrees
     * of freedom per face over all elements
     * of this collection.
     */
    unsigned int max_dofs_per_face () const;

    /**
     * Return the maximal number of degrees
     * of freedom per cell over all elements
     * of this collection.
     */
    unsigned int max_dofs_per_cell () const;

    /**
     * Return an estimate for the memory
     * allocated for this object.
     */
    std::size_t memory_consumption () const;


    /**
     * Return whether all elements
     * in this collection
     * implement the hanging node
     * constraints in the new way,
     * which has to be used to make
     * elements "hp compatible".
     * If this is not the case,
     * the function returns false,
     * which implies, that at least
     * one element in the FECollection
     * does not support the new face
     * interface constraints.
     * On the other hand, if this
     * method does return
     * true, this does not imply
     * that the hp method will work!
     *
     * This behaviour is related to
     * the fact, that FiniteElement
     * classes, which provide the
     * new style hanging node constraints
     * might still not provide
     * them for all possible cases.
     * If FE_Q and FE_RaviartThomas
     * elements are included in the
     * FECollection and both properly implement
     * the get_face_interpolation_matrix
     * method, this method will return
     * true. But the get_face_interpolation_matrix
     * might still fail to find an interpolation
     * matrix between these two elements.
     */
    bool hp_constraints_are_implemented () const;

    /**
     * Return a component mask with as many elements as this
     * object has vector components and of which exactly the
     * one component is true that corresponds to the given
     * argument.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments.
     * It verifies that it gets the same result from every one
     * of the elements that are stored in this FECollection. If
     * this is not the case, it throws an exception.
     *
     * @param scalar An object that represents a single scalar
     * vector component of this finite element.
     * @return A component mask that is false in all components
     * except for the one that corresponds to the argument.
     */
    ComponentMask
    component_mask (const FEValuesExtractors::Scalar &scalar) const;

    /**
     * Return a component mask with as many elements as this
     * object has vector components and of which exactly the
     * <code>dim</code> components are true that correspond to the given
     * argument.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments.
     * It verifies that it gets the same result from every one
     * of the elements that are stored in this FECollection. If
     * this is not the case, it throws an exception.
     *
     * @param vector An object that represents dim
     * vector components of this finite element.
     * @return A component mask that is false in all components
     * except for the ones that corresponds to the argument.
     */
    ComponentMask
    component_mask (const FEValuesExtractors::Vector &vector) const;

    /**
     * Return a component mask with as many elements as this
     * object has vector components and of which exactly the
     * <code>dim*(dim+1)/2</code> components are true that
     * correspond to the given argument.
     *
     * @note This function is the equivalent of
     * FiniteElement::component_mask() with the same arguments.
     * It verifies that it gets the same result from every one
     * of the elements that are stored in this FECollection. If
     * this is not the case, it throws an exception.
     *
     * @param sym_tensor An object that represents dim*(dim+1)/2
     * components of this finite element that are jointly to be
     * interpreted as forming a symmetric tensor.
     * @return A component mask that is false in all components
     * except for the ones that corresponds to the argument.
     */
    ComponentMask
    component_mask (const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

    /**
    * Given a block mask (see @ref GlossBlockMask "this glossary entry"),
    * produce a component mask (see @ref GlossComponentMask "this glossary entry")
    * that represents the components that correspond to the blocks selected in
    * the input argument. This is essentially a conversion operator from
    * BlockMask to ComponentMask.
    *
    * @note This function is the equivalent of
    * FiniteElement::component_mask() with the same arguments.
    * It verifies that it gets the same result from every one
    * of the elements that are stored in this FECollection. If
    * this is not the case, it throws an exception.
    *
    * @param block_mask The mask that selects individual blocks of the finite
    * element
    * @return A mask that selects those components corresponding to the selected
    * blocks of the input argument.
    */
    ComponentMask
    component_mask (const BlockMask &block_mask) const;

    /**
    * Return a block mask with as many elements as this
    * object has blocks and of which exactly the
    * one component is true that corresponds to the given
    * argument. See @ref GlossBlockMask "the glossary"
    * for more information.
    *
    * @note This function will only succeed if the scalar referenced
    * by the argument encompasses a complete block. In other words,
    * if, for example, you pass an extractor for the single
    * $x$ velocity and this object represents an FE_RaviartThomas
    * object, then the single scalar object you selected is part
    * of a larger block and consequently there is no block mask that
    * would represent it. The function will then produce an exception.
    *
    * @note This function is the equivalent of
    * FiniteElement::component_mask() with the same arguments.
    * It verifies that it gets the same result from every one
    * of the elements that are stored in this FECollection. If
    * this is not the case, it throws an exception.
    *
    * @param scalar An object that represents a single scalar
    * vector component of this finite element.
    * @return A component mask that is false in all components
    * except for the one that corresponds to the argument.
    */
    BlockMask
    block_mask (const FEValuesExtractors::Scalar &scalar) const;

    /**
    * Return a component mask with as many elements as this
    * object has vector components and of which exactly the
    * <code>dim</code> components are true that correspond to the given
    * argument. See @ref GlossBlockMask "the glossary"
    * for more information.
    *
    * @note This function is the equivalent of
    * FiniteElement::component_mask() with the same arguments.
    * It verifies that it gets the same result from every one
    * of the elements that are stored in this FECollection. If
    * this is not the case, it throws an exception.
    *
    * @note The same caveat applies as to the version of the function above:
    * The extractor object passed as argument must be so that it corresponds
    * to full blocks and does not split blocks of this element.
    *
    * @param vector An object that represents dim
    * vector components of this finite element.
    * @return A component mask that is false in all components
    * except for the ones that corresponds to the argument.
    */
    BlockMask
    block_mask (const FEValuesExtractors::Vector &vector) const;

    /**
    * Return a component mask with as many elements as this
    * object has vector components and of which exactly the
    * <code>dim*(dim+1)/2</code> components are true that
    * correspond to the given argument. See @ref GlossBlockMask "the glossary"
    * for more information.
    *
    * @note The same caveat applies as to the version of the function above:
    * The extractor object passed as argument must be so that it corresponds
    * to full blocks and does not split blocks of this element.
    *
    * @note This function is the equivalent of
    * FiniteElement::component_mask() with the same arguments.
    * It verifies that it gets the same result from every one
    * of the elements that are stored in this FECollection. If
    * this is not the case, it throws an exception.
    *
    * @param sym_tensor An object that represents dim*(dim+1)/2
    * components of this finite element that are jointly to be
    * interpreted as forming a symmetric tensor.
    * @return A component mask that is false in all components
    * except for the ones that corresponds to the argument.
    */
    BlockMask
    block_mask (const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

    /**
    * Given a component mask (see @ref GlossComponentMask "this glossary entry"),
    * produce a block mask (see @ref GlossBlockMask "this glossary entry")
    * that represents the blocks that correspond to the components selected in
    * the input argument. This is essentially a conversion operator from
    * ComponentMask to BlockMask.
    *
    * @note This function will only succeed if the components referenced
    * by the argument encompasses complete blocks. In other words,
    * if, for example, you pass an component mask for the single
    * $x$ velocity and this object represents an FE_RaviartThomas
    * object, then the single component you selected is part
    * of a larger block and consequently there is no block mask that
    * would represent it. The function will then produce an exception.
    *
    * @note This function is the equivalent of
    * FiniteElement::component_mask() with the same arguments.
    * It verifies that it gets the same result from every one
    * of the elements that are stored in this FECollection. If
    * this is not the case, it throws an exception.
    *
    * @param component_mask The mask that selects individual components of the finite
    * element
    * @return A mask that selects those blocks corresponding to the selected
    * blocks of the input argument.
    */
    BlockMask
    block_mask (const ComponentMask &component_mask) const;


    /**
     * Exception
     */
    DeclException0 (ExcNoFiniteElements);

  private:
    /**
     * Array of pointers to the finite
     * elements stored by this collection.
     */
    std::vector<std_cxx11::shared_ptr<const FiniteElement<dim,spacedim> > > finite_elements;
  };



  /* --------------- inline functions ------------------- */

  template <int dim, int spacedim>
  inline
  unsigned int
  FECollection<dim,spacedim>::size () const
  {
    return finite_elements.size();
  }


  template <int dim, int spacedim>
  inline
  unsigned int
  FECollection<dim,spacedim>::n_components () const
  {
    Assert (finite_elements.size () > 0, ExcNoFiniteElements());

    // note that there is no need
    // here to enforce that indeed
    // all elements have the same
    // number of components since we
    // have already done this when
    // adding a new element to the
    // collection.

    return finite_elements[0]->n_components ();
  }


  template <int dim, int spacedim>
  inline
  const FiniteElement<dim,spacedim> &
  FECollection<dim,spacedim>::operator[] (const unsigned int index) const
  {
    Assert (index < finite_elements.size(),
            ExcIndexRange (index, 0, finite_elements.size()));
    return *finite_elements[index];
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_vertex () const
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_vertex > max)
        max = finite_elements[i]->dofs_per_vertex;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_line () const
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_line > max)
        max = finite_elements[i]->dofs_per_line;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_quad () const
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_quad > max)
        max = finite_elements[i]->dofs_per_quad;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_hex () const
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_hex > max)
        max = finite_elements[i]->dofs_per_hex;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_face () const
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_face > max)
        max = finite_elements[i]->dofs_per_face;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_cell () const
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_cell > max)
        max = finite_elements[i]->dofs_per_cell;

    return max;
  }


  template <int dim, int spacedim>
  bool
  FECollection<dim,spacedim>::hp_constraints_are_implemented () const
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());

    bool hp_constraints = true;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      hp_constraints = hp_constraints &&
                       finite_elements[i]->hp_constraints_are_implemented();

    return hp_constraints;
  }


} // namespace hp

DEAL_II_NAMESPACE_CLOSE

#endif
