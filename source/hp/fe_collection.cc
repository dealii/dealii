// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


#include <deal.II/base/memory_consumption.h>

#include <deal.II/hp/fe_collection.h>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::find_least_face_dominating_fe(
    const std::set<unsigned int> &fes) const
  {
    return find_dominated_fe(find_common_fes(fes, /*codim*/ 1),
                             /*codim*/ 1);
  }



  template <int dim, int spacedim>
  std::set<unsigned int>
  FECollection<dim, spacedim>::find_common_fes(
    const std::set<unsigned int> &fes,
    const unsigned int            codim) const
  {
    // Validate user inputs.
    Assert(codim <= dim, ExcImpossibleInDim(dim));
    for (const auto &fe : fes)
      {
        (void)fe;
        AssertIndexRange(fe, finite_elements.size());
      }

    // Check if any element of this FECollection is able to dominate all
    // elements of @p fes. If one was found, we add it to the set of
    // dominating elements.
    std::set<unsigned int> dominating_fes;
    for (unsigned int current_fe = 0; current_fe < finite_elements.size();
         ++current_fe)
      {
        // Check if current_fe can dominate all elements in @p fes.
        FiniteElementDomination::Domination domination =
          FiniteElementDomination::no_requirements;
        for (const auto &other_fe : fes)
          domination =
            domination & finite_elements[current_fe]->compare_for_domination(
                           *finite_elements[other_fe], codim);

        // If current_fe dominates, add it to the set.
        if ((domination == FiniteElementDomination::this_element_dominates) ||
            (domination == FiniteElementDomination::either_element_can_dominate
             /*covers cases like {Q2,Q3,Q1,Q1} with fes={2,3}*/))
          dominating_fes.insert(current_fe);
      }
    return dominating_fes;
  }



  template <int dim, int spacedim>
  std::set<unsigned int>
  FECollection<dim, spacedim>::find_enclosing_fes(
    const std::set<unsigned int> &fes,
    const unsigned int            codim) const
  {
    // Validate user inputs.
    Assert(codim <= dim, ExcImpossibleInDim(dim));
    for (const auto &fe : fes)
      {
        (void)fe;
        AssertIndexRange(fe, finite_elements.size());
      }

    // Check if any element of this FECollection is dominated by all
    // elements of @p fes. If one was found, we add it to the set of
    // dominated elements.
    std::set<unsigned int> dominated_fes;
    for (unsigned int current_fe = 0; current_fe < finite_elements.size();
         ++current_fe)
      {
        // Check if current_fe is dominated by all other elements in @p fes.
        FiniteElementDomination::Domination domination =
          FiniteElementDomination::no_requirements;
        for (const auto &other_fe : fes)
          domination =
            domination & finite_elements[current_fe]->compare_for_domination(
                           *finite_elements[other_fe], codim);

        // If current_fe is dominated, add it to the set.
        if ((domination == FiniteElementDomination::other_element_dominates) ||
            (domination == FiniteElementDomination::either_element_can_dominate
             /*covers cases like {Q2,Q3,Q1,Q1} with fes={2,3}*/))
          dominated_fes.insert(current_fe);
      }
    return dominated_fes;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::find_dominating_fe(
    const std::set<unsigned int> &fes,
    const unsigned int            codim) const
  {
    // Validate user inputs.
    Assert(codim <= dim, ExcImpossibleInDim(dim));
    for (const auto &fe : fes)
      {
        (void)fe;
        AssertIndexRange(fe, finite_elements.size());
      }

    // If the set of elements contains only a single element,
    // then this very element is considered to be the dominating one.
    if (fes.size() == 1)
      return *fes.begin();

    // There may also be others, in which case we'll check if any of these
    // elements is able to dominate all others. If one was found, we stop
    // looking further and return the dominating element.
    for (const auto &current_fe : fes)
      {
        // Check if current_fe can dominate all elements in @p fes.
        FiniteElementDomination::Domination domination =
          FiniteElementDomination::no_requirements;
        for (const auto &other_fe : fes)
          if (current_fe != other_fe)
            domination =
              domination & finite_elements[current_fe]->compare_for_domination(
                             *finite_elements[other_fe], codim);

        // If current_fe dominates, return its index.
        if ((domination == FiniteElementDomination::this_element_dominates) ||
            (domination == FiniteElementDomination::either_element_can_dominate
             /*covers cases like {Q2,Q3,Q1,Q1} with fes={2,3}*/))
          return current_fe;
      }

    // If we couldn't find the dominating object, return an invalid one.
    return numbers::invalid_unsigned_int;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::find_dominated_fe(
    const std::set<unsigned int> &fes,
    const unsigned int            codim) const
  {
    // Validate user inputs.
    Assert(codim <= dim, ExcImpossibleInDim(dim));
    for (const auto &fe : fes)
      {
        (void)fe;
        AssertIndexRange(fe, finite_elements.size());
      }

    // If the set of elements contains only a single element,
    // then this very element is considered to be the dominated one.
    if (fes.size() == 1)
      return *fes.begin();

    // There may also be others, in which case we'll check if any of these
    // elements is dominated by all others. If one was found, we stop
    // looking further and return the dominated element.
    for (const auto &current_fe : fes)
      {
        // Check if current_fe is dominated by all other elements in @p fes.
        FiniteElementDomination::Domination domination =
          FiniteElementDomination::no_requirements;
        for (const auto &other_fe : fes)
          if (current_fe != other_fe)
            domination =
              domination & finite_elements[current_fe]->compare_for_domination(
                             *finite_elements[other_fe], codim);

        // If current_fe is dominated, return its index.
        if ((domination == FiniteElementDomination::other_element_dominates) ||
            (domination == FiniteElementDomination::either_element_can_dominate
             /*covers cases like {Q2,Q3,Q1,Q1} with fes={2,3}*/))
          return current_fe;
      }

    // If we couldn't find the dominated object, return an invalid one.
    return numbers::invalid_unsigned_int;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::find_dominating_fe_extended(
    const std::set<unsigned int> &fes,
    const unsigned int            codim) const
  {
    unsigned int fe_index = find_dominating_fe(fes, codim);

    if (fe_index == numbers::invalid_unsigned_int)
      {
        const std::set<unsigned int> dominating_fes =
          find_common_fes(fes, codim);
        fe_index = find_dominated_fe(dominating_fes, codim);
      }

    return fe_index;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::find_dominated_fe_extended(
    const std::set<unsigned int> &fes,
    const unsigned int            codim) const
  {
    unsigned int fe_index = find_dominated_fe(fes, codim);

    if (fe_index == numbers::invalid_unsigned_int)
      {
        const std::set<unsigned int> dominated_fes =
          find_enclosing_fes(fes, codim);
        fe_index = find_dominating_fe(dominated_fes, codim);
      }

    return fe_index;
  }



  template <int dim, int spacedim>
  FECollection<dim, spacedim>::FECollection()
  {
    set_default_hierarchy();
  }



  template <int dim, int spacedim>
  FECollection<dim, spacedim>::FECollection(
    const FiniteElement<dim, spacedim> &fe)
    : FECollection()
  {
    push_back(fe);
  }



  template <int dim, int spacedim>
  FECollection<dim, spacedim>::FECollection(
    const std::vector<const FiniteElement<dim, spacedim> *> &fes)
    : FECollection()
  {
    Assert(fes.size() > 0,
           ExcMessage("Need to pass at least one finite element."));

    for (unsigned int i = 0; i < fes.size(); ++i)
      push_back(*fes[i]);
  }



  template <int dim, int spacedim>
  void
  FECollection<dim, spacedim>::push_back(
    const FiniteElement<dim, spacedim> &new_fe)
  {
    // check that the new element has the right
    // number of components. only check with
    // the first element, since all the other
    // elements have already passed the test
    // against the first element
    if (finite_elements.size() != 0)
      Assert(new_fe.n_components() == finite_elements[0]->n_components(),
             ExcMessage("All elements inside a collection need to have the "
                        "same number of vector components!"));

    finite_elements.push_back(new_fe.clone());
  }



  template <int dim, int spacedim>
  void
  FECollection<dim, spacedim>::set_hierarchy(
    const std::function<
      unsigned int(const typename hp::FECollection<dim, spacedim> &,
                   const unsigned int)> &next,
    const std::function<
      unsigned int(const typename hp::FECollection<dim, spacedim> &,
                   const unsigned int)> &prev)
  {
    // copy hierarchy functions
    hierarchy_next = next;
    hierarchy_prev = prev;
  }



  template <int dim, int spacedim>
  void
  FECollection<dim, spacedim>::set_default_hierarchy()
  {
    // establish hierarchy corresponding to order of indices
    set_hierarchy(&DefaultHierarchy::next_index,
                  &DefaultHierarchy::previous_index);
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::next_in_hierarchy(
    const unsigned int fe_index) const
  {
    AssertIndexRange(fe_index, size());

    const unsigned int new_fe_index = hierarchy_next(*this, fe_index);
    AssertIndexRange(new_fe_index, size());

    return new_fe_index;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::previous_in_hierarchy(
    const unsigned int fe_index) const
  {
    AssertIndexRange(fe_index, size());

    const unsigned int new_fe_index = hierarchy_prev(*this, fe_index);
    AssertIndexRange(new_fe_index, size());

    return new_fe_index;
  }



  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(
    const FEValuesExtractors::Scalar &scalar) const
  {
    Assert(size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(scalar);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      Assert(mask == (*this)[c].component_mask(scalar), ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(
    const FEValuesExtractors::Vector &vector) const
  {
    Assert(size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(vector);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      Assert(mask == (*this)[c].component_mask(vector), ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(
    const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
  {
    Assert(size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(sym_tensor);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      Assert(mask == (*this)[c].component_mask(sym_tensor), ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(const BlockMask &block_mask) const
  {
    Assert(size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(block_mask);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      Assert(mask == (*this)[c].component_mask(block_mask),
             ExcMessage("Not all elements of this collection agree on what "
                        "the appropriate mask should be."));

    return mask;
  }


  template <int dim, int spacedim>
  BlockMask
  FECollection<dim, spacedim>::block_mask(
    const FEValuesExtractors::Scalar &scalar) const
  {
    Assert(size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(scalar);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      Assert(mask == (*this)[c].block_mask(scalar),
             ExcMessage("Not all elements of this collection agree on what "
                        "the appropriate mask should be."));

    return mask;
  }


  template <int dim, int spacedim>
  BlockMask
  FECollection<dim, spacedim>::block_mask(
    const FEValuesExtractors::Vector &vector) const
  {
    Assert(size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(vector);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      Assert(mask == (*this)[c].block_mask(vector),
             ExcMessage("Not all elements of this collection agree on what "
                        "the appropriate mask should be."));

    return mask;
  }


  template <int dim, int spacedim>
  BlockMask
  FECollection<dim, spacedim>::block_mask(
    const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
  {
    Assert(size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(sym_tensor);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      Assert(mask == (*this)[c].block_mask(sym_tensor),
             ExcMessage("Not all elements of this collection agree on what "
                        "the appropriate mask should be."));

    return mask;
  }



  template <int dim, int spacedim>
  BlockMask
  FECollection<dim, spacedim>::block_mask(
    const ComponentMask &component_mask) const
  {
    Assert(size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(component_mask);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      Assert(mask == (*this)[c].block_mask(component_mask),
             ExcMessage("Not all elements of this collection agree on what "
                        "the appropriate mask should be."));

    return mask;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::n_blocks() const
  {
    Assert(finite_elements.size() > 0, ExcNoFiniteElements());

    const unsigned int nb = finite_elements[0]->n_blocks();
    for (unsigned int i = 1; i < finite_elements.size(); ++i)
      Assert(finite_elements[i]->n_blocks() == nb,
             ExcMessage("Not all finite elements in this collection have "
                        "the same number of components."));

    return nb;
  }



  template <int dim, int spacedim>
  std::size_t
  FECollection<dim, spacedim>::memory_consumption() const
  {
    std::size_t mem =
      (sizeof(*this) + MemoryConsumption::memory_consumption(finite_elements));
    for (unsigned int i = 0; i < finite_elements.size(); ++i)
      mem += finite_elements[i]->memory_consumption();

    return mem;
  }
} // namespace hp



// explicit instantiations
#include "fe_collection.inst"


DEAL_II_NAMESPACE_CLOSE
