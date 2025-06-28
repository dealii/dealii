// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/memory_consumption.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>

#include <deque>
#include <limits>
#include <set>



DEAL_II_NAMESPACE_OPEN

namespace hp
{
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
    // check that the new element has the right number of components. only check
    // with the first element, since all the other elements have already passed
    // the test against the first element
    Assert(this->empty() ||
             new_fe.n_components() == this->operator[](0).n_components(),
           ExcMessage("All elements inside a collection need to have the "
                      "same number of vector components!"));

    Collection<FiniteElement<dim, spacedim>>::push_back(new_fe.clone());
  }



  template <int dim, int spacedim>
  const MappingCollection<dim, spacedim> &
  FECollection<dim, spacedim>::get_reference_cell_default_linear_mapping() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    // Since we can only add elements to an FECollection, we are safe comparing
    // the sizes of this object and the MappingCollection. One circumstance that
    // might lead to their sizes diverging is this:
    // - An FECollection is created and then this function is called. The shared
    //   map is now initialized.
    // - A second FECollection is made as a copy of this one. The shared map is
    //   not recreated.
    // - The second FECollection is then resized by adding a new FE. The shared
    //   map is thus invalid for the second instance.
    if (!reference_cell_default_linear_mapping ||
        reference_cell_default_linear_mapping->size() != this->size())
      {
        auto &this_nc = const_cast<FECollection<dim, spacedim> &>(*this);

        this_nc.reference_cell_default_linear_mapping =
          std::make_shared<MappingCollection<dim, spacedim>>();

        for (const auto &fe : *this)
          this_nc.reference_cell_default_linear_mapping->push_back(
            fe.reference_cell()
              .template get_default_linear_mapping<dim, spacedim>());
      }

    return *reference_cell_default_linear_mapping;
  }



  template <int dim, int spacedim>
  std::set<unsigned int>
  FECollection<dim, spacedim>::find_common_fes(
    const std::set<unsigned int> &fes,
    const unsigned int            codim) const
  {
    if constexpr (running_in_debug_mode())
      {
        // Validate user inputs.
        Assert(codim <= dim, ExcImpossibleInDim(dim));
        Assert(this->size() > 0, ExcEmptyObject());
        for (const auto &fe : fes)
          AssertIndexRange(fe, this->size());
      }

    // Check if any element of this FECollection is able to dominate all
    // elements of @p fes. If one was found, we add it to the set of
    // dominating elements.
    std::set<unsigned int> dominating_fes;
    for (unsigned int current_fe = 0; current_fe < this->size(); ++current_fe)
      {
        // Check if current_fe can dominate all elements in @p fes.
        FiniteElementDomination::Domination domination =
          FiniteElementDomination::no_requirements;
        for (const auto &other_fe : fes)
          domination =
            domination &
            this->operator[](current_fe)
              .compare_for_domination(this->operator[](other_fe), codim);

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
    if constexpr (running_in_debug_mode())
      {
        // Validate user inputs.
        Assert(codim <= dim, ExcImpossibleInDim(dim));
        Assert(this->size() > 0, ExcEmptyObject());
        for (const auto &fe : fes)
          AssertIndexRange(fe, this->size());
      }

    // Check if any element of this FECollection is dominated by all
    // elements of @p fes. If one was found, we add it to the set of
    // dominated elements.
    std::set<unsigned int> dominated_fes;
    for (unsigned int current_fe = 0; current_fe < this->size(); ++current_fe)
      {
        // Check if current_fe is dominated by all other elements in @p fes.
        FiniteElementDomination::Domination domination =
          FiniteElementDomination::no_requirements;
        for (const auto &other_fe : fes)
          domination =
            domination &
            this->operator[](current_fe)
              .compare_for_domination(this->operator[](other_fe), codim);

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
    // If the set of elements contains only a single element,
    // then this very element is considered to be the dominating one.
    if (fes.size() == 1)
      return *fes.begin();

    if constexpr (running_in_debug_mode())
      {
        // Validate user inputs.
        Assert(codim <= dim, ExcImpossibleInDim(dim));
        Assert(this->size() > 0, ExcEmptyObject());
        for (const auto &fe : fes)
          AssertIndexRange(fe, this->size());
      }

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
              domination &
              this->operator[](current_fe)
                .compare_for_domination(this->operator[](other_fe), codim);

        // If current_fe dominates, return its index.
        if ((domination == FiniteElementDomination::this_element_dominates) ||
            (domination == FiniteElementDomination::either_element_can_dominate
             /*covers cases like {Q2,Q3,Q1,Q1} with fes={2,3}*/))
          return current_fe;
      }

    // If we couldn't find the dominating object, return an invalid one.
    return numbers::invalid_fe_index;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::find_dominated_fe(
    const std::set<unsigned int> &fes,
    const unsigned int            codim) const
  {
    // If the set of elements contains only a single element,
    // then this very element is considered to be the dominated one.
    if (fes.size() == 1)
      return *fes.begin();

    if constexpr (running_in_debug_mode())
      {
        // Validate user inputs.
        Assert(codim <= dim, ExcImpossibleInDim(dim));
        Assert(this->size() > 0, ExcEmptyObject());
        for (const auto &fe : fes)
          AssertIndexRange(fe, this->size());
      }

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
              domination &
              this->operator[](current_fe)
                .compare_for_domination(this->operator[](other_fe), codim);

        // If current_fe is dominated, return its index.
        if ((domination == FiniteElementDomination::other_element_dominates) ||
            (domination == FiniteElementDomination::either_element_can_dominate
             /*covers cases like {Q2,Q3,Q1,Q1} with fes={2,3}*/))
          return current_fe;
      }

    // If we couldn't find the dominated object, return an invalid one.
    return numbers::invalid_fe_index;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::find_dominating_fe_extended(
    const std::set<unsigned int> &fes,
    const unsigned int            codim) const
  {
    unsigned int fe_index = find_dominating_fe(fes, codim);

    if (fe_index == numbers::invalid_fe_index)
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

    if (fe_index == numbers::invalid_fe_index)
      {
        const std::set<unsigned int> dominated_fes =
          find_enclosing_fes(fes, codim);
        fe_index = find_dominating_fe(dominated_fes, codim);
      }

    return fe_index;
  }



  namespace
  {
    /**
     * Implement the action of the hp_*_dof_identities() functions
     * in a generic way.
     */
    std::vector<std::map<unsigned int, unsigned int>>
    compute_hp_dof_identities(
      const std::set<unsigned int> &fes,
      const std::function<std::vector<std::pair<unsigned int, unsigned int>>(
        const unsigned int,
        const unsigned int)>       &query_identities)
    {
      // Let's deal with the easy cases first. If the set of fe indices is empty
      // or has only one entry, then there are no identities:
      if (fes.size() <= 1)
        return {};

      // If the set has two entries, then the
      // FiniteElement::hp_*_dof_identities() function directly returns what we
      // need. We just need to prefix its output with the respective fe indices:
      if (fes.size() == 2)
        {
          const unsigned int fe_index_1 = *fes.begin();
          const unsigned int fe_index_2 = *(++fes.begin());
          const auto         reduced_identities =
            query_identities(fe_index_1, fe_index_2);

          std::vector<std::map<unsigned int, unsigned int>> complete_identities;

          for (const auto &reduced_identity : reduced_identities)
            {
              // Each identity returned by query_identities() is a pair of
              // dof indices. Prefix each with its fe index and put the result
              // into a vector
              std::map<unsigned int, unsigned int> complete_identity = {
                {fe_index_1, reduced_identity.first},
                {fe_index_2, reduced_identity.second}};
              complete_identities.emplace_back(std::move(complete_identity));
            }

          return complete_identities;
        }

      // Now for the general case of three or more elements:
      //
      // Consider all degrees of freedom of the identified elements (represented
      // via (fe_index,dof_index) pairs) as the nodes in a graph. Then draw
      // edges for all DoFs that are identified based on what the elements
      // selected in the argument say. Let us first build this graph, where we
      // only store the edges of the graph, and as a consequence ignore nodes
      // (DoFs) that simply don't show up at all in any of the identities:
      using Node  = std::pair<unsigned int, unsigned int>;
      using Edge  = std::pair<Node, Node>;
      using Graph = std::set<Edge>;

      Graph identities_graph;
      for (const unsigned int fe_index_1 : fes)
        for (const unsigned int fe_index_2 : fes)
          if (fe_index_1 != fe_index_2)
            for (const auto &identity :
                 query_identities(fe_index_1, fe_index_2))
              identities_graph.emplace(Node(fe_index_1, identity.first),
                                       Node(fe_index_2, identity.second));

      if constexpr (running_in_debug_mode())
        {
          // Now verify that indeed the graph is symmetric: If one element
          // declares that certain ones of its DoFs are to be unified with those
          // of the other, then the other one should agree with this. As a
          // consequence of this test succeeding, we know that the graph is
          // actually undirected.
          for (const auto &edge : identities_graph)
            Assert(identities_graph.find({edge.second, edge.first}) !=
                     identities_graph.end(),
                   ExcInternalError());
        }

      // The next step is that we ought to verify that if there is an identity
      // between (fe1,dof1) and (fe2,dof2), as well as with (fe2,dof2) and
      // (fe3,dof3), then there should also be an identity between (fe1,dof1)
      // and (fe3,dof3). The same logic should apply to chains of four
      // identities.
      //
      // This means that the graph we have built above must be composed of a
      // collection of complete sub-graphs (complete = each possible edge in the
      // sub-graph exists) -- or, using a different term, that the graph
      // consists of a number of "cliques". Each of these cliques is then one
      // extended identity between two or more DoFs, and these are the ones that
      // we will want to return.
      //
      // To ascertain that this is true, and to figure out what we want to
      // return, we need to divide the graph into its sub-graphs and then ensure
      // that each sub-graph is indeed a clique. This is slightly cumbersome,
      // but can be done as follows:
      // - pick one edge 'e' of G
      // - add e=(n1,n2) to the sub-graph SG
      // - set N={n1,n2}
      // - loop over the remainder of the edges 'e' of the graph:
      //   - if 'e' has one or both nodes in N:
      //     - add 'e' to SG and
      //     - add its two nodes to N (they may already be in there)
      //     - remove 'e' from G
      //
      // In general, this may not find the entire sub-graph if the edges are
      // stored in random order. For example, if the graph consisted of the
      // following edges in this order:
      //   (a,b)
      //   (c,d)
      //   (a,c)
      //   (a,d)
      //   (b,c)
      //   (b,d)
      // then the graph itself is a clique, but the algorithm outlined above
      // would skip the edge (c,d) because neither of the nodes are already
      // in the set N which at that point is still (a,b).
      //
      // But, we store the graph with a std::set, which stored edges in sorted
      // order where the order is the lexicographic order of nodes. This ensures
      // that we really capture all edges that correspond to a sub-graph (but
      // we will assert this as well).
      //
      // (For programmatic reasons, we skip the removal of 'e' from G in a first
      // run through because it modifies the graph and thus invalidates
      // iterators. But because SG stores all of these edges, we can remove them
      // all from G after collecting the edges in SG.)
      std::vector<std::map<unsigned int, unsigned int>> identities;
      while (identities_graph.size() > 0)
        {
          Graph          sub_graph;       // SG
          std::set<Node> sub_graph_nodes; // N

          sub_graph.emplace(*identities_graph.begin());
          sub_graph_nodes.emplace(identities_graph.begin()->first);
          sub_graph_nodes.emplace(identities_graph.begin()->second);

          for (const Edge &e : identities_graph)
            if ((sub_graph_nodes.find(e.first) != sub_graph_nodes.end()) ||
                (sub_graph_nodes.find(e.second) != sub_graph_nodes.end()))
              {
                sub_graph.insert(e);
                sub_graph_nodes.insert(e.first);
                sub_graph_nodes.insert(e.second);
              }

          // We have now obtained a sub-graph from the overall graph.
          // Now remove it from the bigger graph
          for (const Edge &e : sub_graph)
            identities_graph.erase(e);

          if constexpr (running_in_debug_mode())
            {
              // There are three checks we ought to perform:
              // - That the sub-graph is undirected, i.e. that every edge
              // appears
              //   in both directions
              for (const auto &edge : sub_graph)
                Assert(sub_graph.find({edge.second, edge.first}) !=
                         sub_graph.end(),
                       ExcInternalError());

              // - None of the nodes in the sub-graph should have appeared in
              //   any of the other sub-graphs. If they did, then we have a bug
              //   in extracting sub-graphs. This is actually more easily
              //   checked the other way around: none of the nodes of the
              //   sub-graph we just extracted should be in any of the edges of
              //   the *remaining* graph
              for (const Node &n : sub_graph_nodes)
                for (const Edge &e : identities_graph)
                  Assert((n != e.first) && (n != e.second), ExcInternalError());
              // - Second, the sub-graph we just extracted needs to be complete,
              //   i.e.,
              //   be a "clique". We check this by counting how many edges it
              //   has. for 'n' nodes in 'N', we need to have n*(n-1) edges (we
              //   store both directed edges).
              Assert(sub_graph.size() ==
                       sub_graph_nodes.size() * (sub_graph_nodes.size() - 1),
                     ExcInternalError());
            }

          // At this point we're sure that we have extracted a complete
          // sub-graph ("clique"). The DoFs involved are all identical then, and
          // we will store this identity so we can return it later.
          //
          // The sub-graph is given as a set of Node objects, which is just
          // a collection of (fe_index,dof_index) pairs. Because each
          // fe_index can only appear once there, a better data structure
          // is a std::map from fe_index to dof_index, which can conveniently
          // be initialized from a range of iterators to pairs:
          identities.emplace_back(sub_graph_nodes.begin(),
                                  sub_graph_nodes.end());
          Assert(identities.back().size() == sub_graph_nodes.size(),
                 ExcInternalError());
        }

      return identities;
    }
  } // namespace



  template <int dim, int spacedim>
  std::vector<std::map<unsigned int, unsigned int>>
  FECollection<dim, spacedim>::hp_vertex_dof_identities(
    const std::set<unsigned int> &fes) const
  {
    auto query_vertex_dof_identities = [this](const unsigned int fe_index_1,
                                              const unsigned int fe_index_2) {
      return (*this)[fe_index_1].hp_vertex_dof_identities((*this)[fe_index_2]);
    };
    return compute_hp_dof_identities(fes, query_vertex_dof_identities);
  }



  template <int dim, int spacedim>
  std::vector<std::map<unsigned int, unsigned int>>
  FECollection<dim, spacedim>::hp_line_dof_identities(
    const std::set<unsigned int> &fes) const
  {
    auto query_line_dof_identities = [this](const unsigned int fe_index_1,
                                            const unsigned int fe_index_2) {
      return (*this)[fe_index_1].hp_line_dof_identities((*this)[fe_index_2]);
    };
    return compute_hp_dof_identities(fes, query_line_dof_identities);
  }



  template <int dim, int spacedim>
  std::vector<std::map<unsigned int, unsigned int>>
  FECollection<dim, spacedim>::hp_quad_dof_identities(
    const std::set<unsigned int> &fes,
    const unsigned int            face_no) const
  {
    auto query_quad_dof_identities = [this,
                                      face_no](const unsigned int fe_index_1,
                                               const unsigned int fe_index_2) {
      return (*this)[fe_index_1].hp_quad_dof_identities((*this)[fe_index_2],
                                                        face_no);
    };
    return compute_hp_dof_identities(fes, query_quad_dof_identities);
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
  std::vector<unsigned int>
  FECollection<dim, spacedim>::get_hierarchy_sequence(
    const unsigned int fe_index) const
  {
    AssertIndexRange(fe_index, this->size());

    std::deque<unsigned int> sequence = {fe_index};

    // get predecessors
    {
      unsigned int front = sequence.front();
      unsigned int previous;
      while ((previous = previous_in_hierarchy(front)) != front)
        {
          sequence.push_front(previous);
          front = previous;

          Assert(sequence.size() <= this->size(),
                 ExcMessage(
                   "The registered hierarchy is not terminated: "
                   "previous_in_hierarchy() does not stop at a final index."));
        }
    }

    // get successors
    {
      unsigned int back = sequence.back();
      unsigned int next;
      while ((next = next_in_hierarchy(back)) != back)
        {
          sequence.push_back(next);
          back = next;

          Assert(sequence.size() <= this->size(),
                 ExcMessage(
                   "The registered hierarchy is not terminated: "
                   "next_in_hierarchy() does not stop at a final index."));
        }
    }

    return {sequence.begin(), sequence.end()};
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::next_in_hierarchy(
    const unsigned int fe_index) const
  {
    AssertIndexRange(fe_index, this->size());

    const unsigned int new_fe_index = hierarchy_next(*this, fe_index);
    AssertIndexRange(new_fe_index, this->size());

    return new_fe_index;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::previous_in_hierarchy(
    const unsigned int fe_index) const
  {
    AssertIndexRange(fe_index, this->size());

    const unsigned int new_fe_index = hierarchy_prev(*this, fe_index);
    AssertIndexRange(new_fe_index, this->size());

    return new_fe_index;
  }



  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(
    const FEValuesExtractors::Scalar &scalar) const
  {
    Assert(this->size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(scalar);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < this->size(); ++c)
      Assert(mask == (*this)[c].component_mask(scalar), ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(
    const FEValuesExtractors::Vector &vector) const
  {
    Assert(this->size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(vector);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < this->size(); ++c)
      Assert(mask == (*this)[c].component_mask(vector), ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(
    const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
  {
    Assert(this->size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(sym_tensor);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < this->size(); ++c)
      Assert(mask == (*this)[c].component_mask(sym_tensor), ExcInternalError());

    return mask;
  }


  template <int dim, int spacedim>
  ComponentMask
  FECollection<dim, spacedim>::component_mask(const BlockMask &block_mask) const
  {
    Assert(this->size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const ComponentMask mask = (*this)[0].component_mask(block_mask);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < this->size(); ++c)
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
    Assert(this->size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(scalar);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < this->size(); ++c)
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
    Assert(this->size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(vector);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < this->size(); ++c)
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
    Assert(this->size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(sym_tensor);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < this->size(); ++c)
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
    Assert(this->size() > 0,
           ExcMessage("This collection contains no finite element."));

    // get the mask from the first element of the collection
    const BlockMask mask = (*this)[0].block_mask(component_mask);

    // but then also verify that the other elements of the collection
    // would return the same mask
    for (unsigned int c = 1; c < this->size(); ++c)
      Assert(mask == (*this)[c].block_mask(component_mask),
             ExcMessage("Not all elements of this collection agree on what "
                        "the appropriate mask should be."));

    return mask;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::n_blocks() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    const unsigned int nb = this->operator[](0).n_blocks();
    for (unsigned int i = 1; i < this->size(); ++i)
      Assert(this->operator[](i).n_blocks() == nb,
             ExcMessage("Not all finite elements in this collection have "
                        "the same number of components."));

    return nb;
  }
} // namespace hp



// explicit instantiations
#include "hp/fe_collection.inst"


DEAL_II_NAMESPACE_CLOSE
