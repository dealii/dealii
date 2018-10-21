// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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
  FECollection<dim, spacedim>::find_least_face_dominating_fe_in_collection(
    const std::set<unsigned int> &fes) const
  {
    for (auto it = fes.cbegin(); it != fes.cend(); ++it)
      DEAL_II_AssertIndexRange(*it, finite_elements.size());

    // If the set of elements to be dominated contains only a single element X,
    // then by definition the dominating set contains this single element X
    // (because each element can dominate itself). There may also be others,
    // say Y1...YN. Next you have to find one or more elements in the dominating
    // set {X,Y1...YN} that is the weakest. Well, you can't find one that is
    // weaker than X because if it were, it would not dominate X. In other
    // words, X is guaranteed to be in the subset of {X,Y1...YN} of weakest
    // dominating elements. Since we only guarantee that the function returns
    // one of them, we may as well return X right away.
    if (fes.size() == 1)
      return *fes.begin();

    std::set<unsigned int> candidate_fes;

    // first loop over all FEs and check which can dominate those given in @p fes:
    for (unsigned int cur_fe = 0; cur_fe < finite_elements.size(); ++cur_fe)
      {
        FiniteElementDomination::Domination domination =
          FiniteElementDomination::no_requirements;
        // check if cur_fe can dominate all FEs in @p fes:
        for (const auto &other_fe : fes)
          domination =
            domination & finite_elements[cur_fe]->compare_for_face_domination(
                           *finite_elements[other_fe]);

        // if we found dominating element, keep them in a set.
        if (
          domination == FiniteElementDomination::this_element_dominates ||
          domination == FiniteElementDomination::either_element_can_dominate /*covers cases like {Q2,Q3,Q1,Q1} with fes={2,3}*/)
          candidate_fes.insert(cur_fe);
      }

    // among the ones we found, pick one that is dominated by all others and
    // thus should represent the largest FE space.
    if (candidate_fes.size() == 1)
      {
        return *candidate_fes.begin();
      }
    else
      for (const auto &current_fe : candidate_fes)
        {
          FiniteElementDomination::Domination domination =
            FiniteElementDomination::no_requirements;

          for (const auto &other_fe : candidate_fes)
            if (current_fe != other_fe)
              domination =
                domination &
                finite_elements[current_fe]->compare_for_face_domination(
                  *finite_elements[other_fe]);

          if ((domination ==
               FiniteElementDomination::other_element_dominates) ||
              (domination ==
               FiniteElementDomination::either_element_can_dominate
               /*covers cases like candidate_fes={Q1,Q1}*/))
            return current_fe;
        }
    // We couldn't find the FE, return invalid_unsigned_int :
    return numbers::invalid_unsigned_int;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::find_least_face_dominating_fe(
    const std::set<unsigned int> &fes) const
  {
    return find_least_face_dominating_fe_in_collection(fes);
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::find_face_dominating_fe_in_subset(
    const std::set<unsigned int> &fes) const
  {
    for (auto it = fes.cbegin(); it != fes.cend(); ++it)
      DEAL_II_AssertIndexRange(*it, finite_elements.size());

    // If the set of elements to be dominated contains only a single element X,
    // then by definition the dominating set contains this single element
    // (because each element can dominate itself).
    // There may also be others, in which case we'll check if any of these
    // elements is able to dominate all others. If one was found, we stop
    // looking further and return the dominating element.
    if (fes.size() == 1)
      return *fes.begin();

    // loop over all finite elements given in the subset
    // and check which one dominates the whole subset
    for (const auto &current_fe : fes)
      {
        FiniteElementDomination::Domination domination =
          FiniteElementDomination::either_element_can_dominate;

        for (const auto &other_fe : fes)
          if (other_fe != current_fe)
            domination =
              domination &
              finite_elements[current_fe]->compare_for_face_domination(
                *finite_elements[other_fe]);

        // see if this element is able to dominate all the other
        // ones, and if so take it
        if ((domination == FiniteElementDomination::this_element_dominates) ||
            (domination ==
             FiniteElementDomination::either_element_can_dominate) ||
            (domination == FiniteElementDomination::no_requirements))
          return current_fe;
      }

    // if we couldn't find the most dominating object
    return numbers::invalid_unsigned_int;
  }



  template <int dim, int spacedim>
  FECollection<dim, spacedim>::FECollection(
    const FiniteElement<dim, spacedim> &fe)
  {
    push_back(fe);
  }



  template <int dim, int spacedim>
  FECollection<dim, spacedim>::FECollection(
    const std::vector<const FiniteElement<dim, spacedim> *> &fes)
  {
    DEAL_II_Assert(fes.size() > 0,
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
      DEAL_II_Assert(new_fe.n_components() ==
                       finite_elements[0]->n_components(),
                     ExcMessage(
                       "All elements inside a collection need to have the "
                       "same number of vector components!"));

    finite_elements.push_back(new_fe.clone());
  }



  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(
    const FEValuesExtractors::Scalar &scalar) const
  {
    DEAL_II_Assert(size() > 0,
                   ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(scalar);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      DEAL_II_Assert(mask == (*this)[c].component_mask(scalar),
                     ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(
    const FEValuesExtractors::Vector &vector) const
  {
    DEAL_II_Assert(size() > 0,
                   ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(vector);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      DEAL_II_Assert(mask == (*this)[c].component_mask(vector),
                     ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(
    const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
  {
    DEAL_II_Assert(size() > 0,
                   ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(sym_tensor);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      DEAL_II_Assert(mask == (*this)[c].component_mask(sym_tensor),
                     ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(const BlockMask &block_mask) const
  {
    DEAL_II_Assert(size() > 0,
                   ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(block_mask);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      DEAL_II_Assert(mask == (*this)[c].component_mask(block_mask),
                     ExcMessage(
                       "Not all elements of this collection agree on what "
                       "the appropriate mask should be."));

    return mask;
  }


  template <int dim, int spacedim>
  BlockMask
  FECollection<dim, spacedim>::block_mask(
    const FEValuesExtractors::Scalar &scalar) const
  {
    DEAL_II_Assert(size() > 0,
                   ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(scalar);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      DEAL_II_Assert(mask == (*this)[c].block_mask(scalar),
                     ExcMessage(
                       "Not all elements of this collection agree on what "
                       "the appropriate mask should be."));

    return mask;
  }


  template <int dim, int spacedim>
  BlockMask
  FECollection<dim, spacedim>::block_mask(
    const FEValuesExtractors::Vector &vector) const
  {
    DEAL_II_Assert(size() > 0,
                   ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(vector);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      DEAL_II_Assert(mask == (*this)[c].block_mask(vector),
                     ExcMessage(
                       "Not all elements of this collection agree on what "
                       "the appropriate mask should be."));

    return mask;
  }


  template <int dim, int spacedim>
  BlockMask
  FECollection<dim, spacedim>::block_mask(
    const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
  {
    DEAL_II_Assert(size() > 0,
                   ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(sym_tensor);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      DEAL_II_Assert(mask == (*this)[c].block_mask(sym_tensor),
                     ExcMessage(
                       "Not all elements of this collection agree on what "
                       "the appropriate mask should be."));

    return mask;
  }



  template <int dim, int spacedim>
  BlockMask
  FECollection<dim, spacedim>::block_mask(
    const ComponentMask &component_mask) const
  {
    DEAL_II_Assert(size() > 0,
                   ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(component_mask);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < size(); ++c)
      DEAL_II_Assert(mask == (*this)[c].block_mask(component_mask),
                     ExcMessage(
                       "Not all elements of this collection agree on what "
                       "the appropriate mask should be."));

    return mask;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::n_blocks() const
  {
    DEAL_II_Assert(finite_elements.size() > 0, ExcNoFiniteElements());

    const unsigned int nb = finite_elements[0]->n_blocks();
    for (unsigned int i = 1; i < finite_elements.size(); ++i)
      DEAL_II_Assert(finite_elements[i]->n_blocks() == nb,
                     ExcMessage(
                       "Not all finite elements in this collection have "
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
