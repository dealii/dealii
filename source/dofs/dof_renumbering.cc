// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q_base.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/king_ordering.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#undef BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace DoFRenumbering
{
  namespace boost
  {
    namespace boosttypes
    {
      using namespace ::boost;

      using Graph     = adjacency_list<vecS,
                                   vecS,
                                   undirectedS,
                                   property<vertex_color_t,
                                            default_color_type,
                                            property<vertex_degree_t, int>>>;
      using Vertex    = graph_traits<Graph>::vertex_descriptor;
      using size_type = graph_traits<Graph>::vertices_size_type;

      using Pair = std::pair<size_type, size_type>;
    } // namespace boosttypes


    namespace internal
    {
      template <int dim, int spacedim>
      void
      create_graph(const DoFHandler<dim, spacedim> &dof_handler,
                   const bool                       use_constraints,
                   boosttypes::Graph               &graph,
                   boosttypes::property_map<boosttypes::Graph,
                                            boosttypes::vertex_degree_t>::type
                     &graph_degree)
      {
        {
          // create intermediate sparsity pattern (faster than directly
          // submitting indices)
          AffineConstraints<double> constraints;
          if (use_constraints)
            DoFTools::make_hanging_node_constraints(dof_handler, constraints);
          constraints.close();
          DynamicSparsityPattern dsp(dof_handler.n_dofs(),
                                     dof_handler.n_dofs());
          DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);

          // submit the entries to the boost graph
          for (unsigned int row = 0; row < dsp.n_rows(); ++row)
            for (unsigned int col = 0; col < dsp.row_length(row); ++col)
              add_edge(row, dsp.column_number(row, col), graph);
        }

        boosttypes::graph_traits<boosttypes::Graph>::vertex_iterator ui, ui_end;

        graph_degree = get(::boost::vertex_degree, graph);
        for (::boost::tie(ui, ui_end) = vertices(graph); ui != ui_end; ++ui)
          graph_degree[*ui] = degree(*ui, graph);
      }
    } // namespace internal


    template <int dim, int spacedim>
    void
    Cuthill_McKee(DoFHandler<dim, spacedim> &dof_handler,
                  const bool                 reversed_numbering,
                  const bool                 use_constraints)
    {
      std::vector<types::global_dof_index> renumbering(
        dof_handler.n_dofs(), numbers::invalid_dof_index);
      compute_Cuthill_McKee(renumbering,
                            dof_handler,
                            reversed_numbering,
                            use_constraints);

      // actually perform renumbering;
      // this is dimension specific and
      // thus needs an own function
      dof_handler.renumber_dofs(renumbering);
    }


    template <int dim, int spacedim>
    void
    compute_Cuthill_McKee(std::vector<types::global_dof_index> &new_dof_indices,
                          const DoFHandler<dim, spacedim>      &dof_handler,
                          const bool reversed_numbering,
                          const bool use_constraints)
    {
      boosttypes::Graph graph(dof_handler.n_dofs());
      boosttypes::property_map<boosttypes::Graph,
                               boosttypes::vertex_degree_t>::type graph_degree;

      internal::create_graph(dof_handler, use_constraints, graph, graph_degree);

      boosttypes::property_map<boosttypes::Graph,
                               boosttypes::vertex_index_t>::type index_map =
        get(::boost::vertex_index, graph);


      std::vector<boosttypes::Vertex> inv_perm(num_vertices(graph));

      if (reversed_numbering == false)
        ::boost::cuthill_mckee_ordering(graph,
                                        inv_perm.rbegin(),
                                        get(::boost::vertex_color, graph),
                                        make_degree_map(graph));
      else
        ::boost::cuthill_mckee_ordering(graph,
                                        inv_perm.begin(),
                                        get(::boost::vertex_color, graph),
                                        make_degree_map(graph));

      for (boosttypes::size_type c = 0; c != inv_perm.size(); ++c)
        new_dof_indices[index_map[inv_perm[c]]] = c;

      Assert(std::find(new_dof_indices.begin(),
                       new_dof_indices.end(),
                       numbers::invalid_dof_index) == new_dof_indices.end(),
             ExcInternalError());
    }



    template <int dim, int spacedim>
    void
    king_ordering(DoFHandler<dim, spacedim> &dof_handler,
                  const bool                 reversed_numbering,
                  const bool                 use_constraints)
    {
      std::vector<types::global_dof_index> renumbering(
        dof_handler.n_dofs(), numbers::invalid_dof_index);
      compute_king_ordering(renumbering,
                            dof_handler,
                            reversed_numbering,
                            use_constraints);

      // actually perform renumbering;
      // this is dimension specific and
      // thus needs an own function
      dof_handler.renumber_dofs(renumbering);
    }


    template <int dim, int spacedim>
    void
    compute_king_ordering(std::vector<types::global_dof_index> &new_dof_indices,
                          const DoFHandler<dim, spacedim>      &dof_handler,
                          const bool reversed_numbering,
                          const bool use_constraints)
    {
      boosttypes::Graph graph(dof_handler.n_dofs());
      boosttypes::property_map<boosttypes::Graph,
                               boosttypes::vertex_degree_t>::type graph_degree;

      internal::create_graph(dof_handler, use_constraints, graph, graph_degree);

      boosttypes::property_map<boosttypes::Graph,
                               boosttypes::vertex_index_t>::type index_map =
        get(::boost::vertex_index, graph);


      std::vector<boosttypes::Vertex> inv_perm(num_vertices(graph));

      if (reversed_numbering == false)
        ::boost::king_ordering(graph, inv_perm.rbegin());
      else
        ::boost::king_ordering(graph, inv_perm.begin());

      for (boosttypes::size_type c = 0; c != inv_perm.size(); ++c)
        new_dof_indices[index_map[inv_perm[c]]] = c;

      Assert(std::find(new_dof_indices.begin(),
                       new_dof_indices.end(),
                       numbers::invalid_dof_index) == new_dof_indices.end(),
             ExcInternalError());
    }



    template <int dim, int spacedim>
    void
    minimum_degree(DoFHandler<dim, spacedim> &dof_handler,
                   const bool                 reversed_numbering,
                   const bool                 use_constraints)
    {
      std::vector<types::global_dof_index> renumbering(
        dof_handler.n_dofs(), numbers::invalid_dof_index);
      compute_minimum_degree(renumbering,
                             dof_handler,
                             reversed_numbering,
                             use_constraints);

      // actually perform renumbering;
      // this is dimension specific and
      // thus needs an own function
      dof_handler.renumber_dofs(renumbering);
    }


    template <int dim, int spacedim>
    void
    compute_minimum_degree(
      std::vector<types::global_dof_index> &new_dof_indices,
      const DoFHandler<dim, spacedim>      &dof_handler,
      const bool                            reversed_numbering,
      const bool                            use_constraints)
    {
      (void)use_constraints;
      Assert(use_constraints == false, ExcNotImplemented());

      // the following code is pretty
      // much a verbatim copy of the
      // sample code for the
      // minimum_degree_ordering manual
      // page from the BOOST Graph
      // Library
      using namespace ::boost;

      int delta = 0;

      // must be BGL directed graph now
      using Graph = adjacency_list<vecS, vecS, directedS>;

      int n = dof_handler.n_dofs();

      Graph G(n);

      std::vector<dealii::types::global_dof_index> dofs_on_this_cell;

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();

          dofs_on_this_cell.resize(dofs_per_cell);

          cell->get_active_or_mg_dof_indices(dofs_on_this_cell);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              if (dofs_on_this_cell[i] > dofs_on_this_cell[j])
                {
                  add_edge(dofs_on_this_cell[i], dofs_on_this_cell[j], G);
                  add_edge(dofs_on_this_cell[j], dofs_on_this_cell[i], G);
                }
        }


      using Vector = std::vector<int>;


      Vector inverse_perm(n, 0);

      Vector perm(n, 0);


      Vector supernode_sizes(n, 1);
      // init has to be 1

      ::boost::property_map<Graph, vertex_index_t>::type id =
        get(vertex_index, G);


      Vector degree(n, 0);


      minimum_degree_ordering(
        G,
        make_iterator_property_map(degree.begin(), id, degree[0]),
        inverse_perm.data(),
        perm.data(),
        make_iterator_property_map(supernode_sizes.begin(),
                                   id,
                                   supernode_sizes[0]),
        delta,
        id);


      for (int i = 0; i < n; ++i)
        {
          Assert(std::find(perm.begin(), perm.end(), i) != perm.end(),
                 ExcInternalError());
          Assert(std::find(inverse_perm.begin(), inverse_perm.end(), i) !=
                   inverse_perm.end(),
                 ExcInternalError());
          Assert(inverse_perm[perm[i]] == i, ExcInternalError());
        }

      if (reversed_numbering == true)
        std::copy(perm.begin(), perm.end(), new_dof_indices.begin());
      else
        std::copy(inverse_perm.begin(),
                  inverse_perm.end(),
                  new_dof_indices.begin());
    }

  } // namespace boost



  template <int dim, int spacedim>
  void
  Cuthill_McKee(DoFHandler<dim, spacedim>                  &dof_handler,
                const bool                                  reversed_numbering,
                const bool                                  use_constraints,
                const std::vector<types::global_dof_index> &starting_indices)
  {
    std::vector<types::global_dof_index> renumbering(
      dof_handler.locally_owned_dofs().n_elements(),
      numbers::invalid_dof_index);
    compute_Cuthill_McKee(renumbering,
                          dof_handler,
                          reversed_numbering,
                          use_constraints,
                          starting_indices);

    // actually perform renumbering;
    // this is dimension specific and
    // thus needs an own function
    dof_handler.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  compute_Cuthill_McKee(
    std::vector<types::global_dof_index>       &new_indices,
    const DoFHandler<dim, spacedim>            &dof_handler,
    const bool                                  reversed_numbering,
    const bool                                  use_constraints,
    const std::vector<types::global_dof_index> &starting_indices,
    const unsigned int                          level)
  {
    const bool reorder_level_dofs =
      (level == numbers::invalid_unsigned_int) ? false : true;

    // see if there is anything to do at all or whether we can skip the work on
    // this processor
    if (dof_handler.locally_owned_dofs().n_elements() == 0)
      {
        Assert(new_indices.empty(), ExcInternalError());
        return;
      }

    // make the connection graph
    //
    // note that if constraints are not requested, then the 'constraints'
    // object will be empty and using it has no effect
    if (reorder_level_dofs == true)
      Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
             ExcDoFHandlerNotInitialized());

    const IndexSet locally_relevant_dofs =
      (reorder_level_dofs == false ?
         DoFTools::extract_locally_relevant_dofs(dof_handler) :
         DoFTools::extract_locally_relevant_level_dofs(dof_handler, level));
    const IndexSet &locally_owned_dofs =
      (reorder_level_dofs == false ? dof_handler.locally_owned_dofs() :
                                     dof_handler.locally_owned_mg_dofs(level));

    AffineConstraints<double> constraints;
    if (use_constraints)
      {
        // reordering with constraints is not yet implemented on a level basis
        Assert(reorder_level_dofs == false, ExcNotImplemented());

        constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
        DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      }
    constraints.close();

    // see if we can get away with the sequential algorithm
    if (locally_owned_dofs.n_elements() == locally_owned_dofs.size())
      {
        AssertDimension(new_indices.size(), locally_owned_dofs.n_elements());

        DynamicSparsityPattern dsp(locally_owned_dofs.size(),
                                   locally_owned_dofs.size());
        if (reorder_level_dofs == false)
          {
            DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
          }
        else
          {
            MGTools::make_sparsity_pattern(dof_handler, dsp, level);
          }

        SparsityTools::reorder_Cuthill_McKee(dsp,
                                             new_indices,
                                             starting_indices);
        if (reversed_numbering)
          new_indices = Utilities::reverse_permutation(new_indices);
      }
    else
      {
        // we are in the parallel case where we need to work in the
        // local index space, i.e., the locally owned part of the
        // sparsity pattern.
        //
        // first figure out whether the user only gave us starting
        // indices that are locally owned, or that are only locally
        // relevant. in the process, also check that all indices
        // really belong to at least the locally relevant ones
        const IndexSet locally_active_dofs =
          (reorder_level_dofs == false ?
             DoFTools::extract_locally_active_dofs(dof_handler) :
             DoFTools::extract_locally_active_level_dofs(dof_handler, level));

        bool needs_locally_active = false;
        for (const auto starting_index : starting_indices)
          {
            if ((needs_locally_active ==
                 /* previously already set to */ true) ||
                (locally_owned_dofs.is_element(starting_index) == false))
              {
                Assert(
                  locally_active_dofs.is_element(starting_index),
                  ExcMessage(
                    "You specified global degree of freedom " +
                    std::to_string(starting_index) +
                    " as a starting index, but this index is not among the "
                    "locally active ones on this processor, as required "
                    "for this function."));
                needs_locally_active = true;
              }
          }

        const IndexSet index_set_to_use =
          (needs_locally_active ? locally_active_dofs : locally_owned_dofs);

        // if this process doesn't own any DoFs (on this level), there is
        // nothing to do
        if (index_set_to_use.n_elements() == 0)
          return;

        // then create first the global sparsity pattern, and then the local
        // sparsity pattern from the global one by transferring its indices to
        // processor-local (locally owned or locally active) index space
        DynamicSparsityPattern dsp(index_set_to_use.size(),
                                   index_set_to_use.size(),
                                   index_set_to_use);
        if (reorder_level_dofs == false)
          {
            DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
          }
        else
          {
            MGTools::make_sparsity_pattern(dof_handler, dsp, level);
          }

        DynamicSparsityPattern local_sparsity(index_set_to_use.n_elements(),
                                              index_set_to_use.n_elements());
        std::vector<types::global_dof_index> row_entries;
        for (unsigned int i = 0; i < index_set_to_use.n_elements(); ++i)
          {
            const types::global_dof_index row =
              index_set_to_use.nth_index_in_set(i);
            const unsigned int row_length = dsp.row_length(row);
            row_entries.clear();
            for (unsigned int j = 0; j < row_length; ++j)
              {
                const unsigned int col = dsp.column_number(row, j);
                if (col != row && index_set_to_use.is_element(col))
                  row_entries.push_back(index_set_to_use.index_within_set(col));
              }
            local_sparsity.add_entries(i,
                                       row_entries.begin(),
                                       row_entries.end(),
                                       true);
          }

        // translate starting indices from global to local indices
        std::vector<types::global_dof_index> local_starting_indices(
          starting_indices.size());
        for (unsigned int i = 0; i < starting_indices.size(); ++i)
          local_starting_indices[i] =
            index_set_to_use.index_within_set(starting_indices[i]);

        // then do the renumbering on the locally owned portion
        AssertDimension(new_indices.size(), locally_owned_dofs.n_elements());
        std::vector<types::global_dof_index> my_new_indices(
          index_set_to_use.n_elements());
        SparsityTools::reorder_Cuthill_McKee(local_sparsity,
                                             my_new_indices,
                                             local_starting_indices);
        if (reversed_numbering)
          my_new_indices = Utilities::reverse_permutation(my_new_indices);

        // now that we have a re-enumeration of all DoFs, we need to throw
        // out the ones that are not locally owned in case we have worked
        // with the locally active ones. that's because the renumbering
        // functions only want new indices for the locally owned DoFs (other
        // processors are responsible for renumbering the ones that are
        // on cell interfaces)
        if (needs_locally_active == true)
          {
            // first step: figure out which DoF indices to eliminate
            IndexSet active_but_not_owned_dofs = locally_active_dofs;
            active_but_not_owned_dofs.subtract_set(locally_owned_dofs);

            std::set<types::global_dof_index> erase_these_indices;
            for (const auto p : active_but_not_owned_dofs)
              {
                const auto index = index_set_to_use.index_within_set(p);
                Assert(index < index_set_to_use.n_elements(),
                       ExcInternalError());
                erase_these_indices.insert(my_new_indices[index]);
                my_new_indices[index] = numbers::invalid_dof_index;
              }
            Assert(erase_these_indices.size() ==
                     active_but_not_owned_dofs.n_elements(),
                   ExcInternalError());
            Assert(static_cast<unsigned int>(
                     std::count(my_new_indices.begin(),
                                my_new_indices.end(),
                                numbers::invalid_dof_index)) ==
                     active_but_not_owned_dofs.n_elements(),
                   ExcInternalError());

            // then compute a renumbering of the remaining ones
            std::vector<types::global_dof_index> translate_indices(
              my_new_indices.size());
            {
              std::set<types::global_dof_index>::const_iterator
                next_erased_index = erase_these_indices.begin();
              types::global_dof_index next_new_index = 0;
              for (unsigned int i = 0; i < translate_indices.size(); ++i)
                if ((next_erased_index != erase_these_indices.end()) &&
                    (*next_erased_index == i))
                  {
                    translate_indices[i] = numbers::invalid_dof_index;
                    ++next_erased_index;
                  }
                else
                  {
                    translate_indices[i] = next_new_index;
                    ++next_new_index;
                  }
              Assert(next_new_index == locally_owned_dofs.n_elements(),
                     ExcInternalError());
            }

            // and then do the renumbering of the result of the
            // Cuthill-McKee algorithm above, right into the output array
            new_indices.clear();
            new_indices.reserve(locally_owned_dofs.n_elements());
            for (const auto &p : my_new_indices)
              if (p != numbers::invalid_dof_index)
                {
                  Assert(translate_indices[p] != numbers::invalid_dof_index,
                         ExcInternalError());
                  new_indices.push_back(translate_indices[p]);
                }
            Assert(new_indices.size() == locally_owned_dofs.n_elements(),
                   ExcInternalError());
          }
        else
          new_indices = std::move(my_new_indices);

        // convert indices back to global index space. in both of the branches
        // above, we ended up with new_indices only containing the local
        // indices of the locally-owned DoFs. so that's where we get the
        // indices
        for (types::global_dof_index &new_index : new_indices)
          new_index = locally_owned_dofs.nth_index_in_set(new_index);
      }
  }



  template <int dim, int spacedim>
  void
  Cuthill_McKee(DoFHandler<dim, spacedim>                  &dof_handler,
                const unsigned int                          level,
                const bool                                  reversed_numbering,
                const std::vector<types::global_dof_index> &starting_indices)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcDoFHandlerNotInitialized());

    std::vector<types::global_dof_index> new_indices(
      dof_handler.locally_owned_mg_dofs(level).n_elements(),
      numbers::invalid_dof_index);

    compute_Cuthill_McKee(new_indices,
                          dof_handler,
                          reversed_numbering,
                          false,
                          starting_indices,
                          level);

    // actually perform renumbering;
    // this is dimension specific and
    // thus needs an own function
    dof_handler.renumber_dofs(level, new_indices);
  }



  template <int dim, int spacedim>
  void
  component_wise(DoFHandler<dim, spacedim>       &dof_handler,
                 const std::vector<unsigned int> &component_order_arg)
  {
    std::vector<types::global_dof_index> renumbering(
      dof_handler.n_locally_owned_dofs(), numbers::invalid_dof_index);

    const types::global_dof_index result =
      compute_component_wise<dim, spacedim>(renumbering,
                                            dof_handler.begin_active(),
                                            dof_handler.end(),
                                            component_order_arg,
                                            false);
    (void)result;

    // If we don't have a renumbering (i.e., when there is 1 component) then
    // return
    if (Utilities::MPI::max(renumbering.size(),
                            dof_handler.get_mpi_communicator()) == 0)
      return;

    // verify that the last numbered
    // degree of freedom is either
    // equal to the number of degrees
    // of freedom in total (the
    // sequential case) or in the
    // distributed case at least
    // makes sense
    Assert((result == dof_handler.n_locally_owned_dofs()) ||
             ((dof_handler.n_locally_owned_dofs() < dof_handler.n_dofs()) &&
              (result <= dof_handler.n_dofs())),
           ExcInternalError());

    dof_handler.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  component_wise(DoFHandler<dim, spacedim>       &dof_handler,
                 const unsigned int               level,
                 const std::vector<unsigned int> &component_order_arg)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcDoFHandlerNotInitialized());

    std::vector<types::global_dof_index> renumbering(
      dof_handler.locally_owned_mg_dofs(level).n_elements(),
      numbers::invalid_dof_index);

    typename DoFHandler<dim, spacedim>::level_cell_iterator start =
      dof_handler.begin(level);
    typename DoFHandler<dim, spacedim>::level_cell_iterator end =
      dof_handler.end(level);

    const types::global_dof_index result =
      compute_component_wise<dim, spacedim>(
        renumbering, start, end, component_order_arg, true);
    (void)result;

    // If we don't have a renumbering (i.e., when there is 1 component) then
    // return
    if (Utilities::MPI::max(renumbering.size(),
                            dof_handler.get_mpi_communicator()) == 0)
      return;

    // verify that the last numbered
    // degree of freedom is either
    // equal to the number of degrees
    // of freedom in total (the
    // sequential case) or in the
    // distributed case at least
    // makes sense
    Assert((result == dof_handler.locally_owned_mg_dofs(level).n_elements()) ||
             ((dof_handler.locally_owned_mg_dofs(level).n_elements() <
               dof_handler.n_dofs(level)) &&
              (result <= dof_handler.n_dofs(level))),
           ExcInternalError());

    dof_handler.renumber_dofs(level, renumbering);
  }



  template <int dim, int spacedim, typename CellIterator>
  types::global_dof_index
  compute_component_wise(std::vector<types::global_dof_index> &new_indices,
                         const CellIterator                   &start,
                         const std_cxx20::type_identity_t<CellIterator> &end,
                         const std::vector<unsigned int> &component_order_arg,
                         const bool                       is_level_operation)
  {
    const hp::FECollection<dim, spacedim> &fe_collection =
      start->get_dof_handler().get_fe_collection();

    // do nothing if the FE has only
    // one component
    if (fe_collection.n_components() == 1)
      {
        new_indices.resize(0);
        return 0;
      }

    // Get a reference to the set of dofs. Note that we assume that all cells
    // are assumed to be on the same level, otherwise the operation doesn't make
    // much sense (we will assert this below).
    const IndexSet &locally_owned_dofs =
      is_level_operation ?
        start->get_dof_handler().locally_owned_mg_dofs(start->level()) :
        start->get_dof_handler().locally_owned_dofs();

    // Copy last argument into a
    // writable vector.
    std::vector<unsigned int> component_order(component_order_arg);
    // If the last argument was an
    // empty vector, set up things to
    // store components in the order
    // found in the system.
    if (component_order.empty())
      for (unsigned int i = 0; i < fe_collection.n_components(); ++i)
        component_order.push_back(i);

    Assert(component_order.size() == fe_collection.n_components(),
           ExcDimensionMismatch(component_order.size(),
                                fe_collection.n_components()));

    for (const unsigned int component : component_order)
      {
        (void)component;
        AssertIndexRange(component, fe_collection.n_components());
      }

    // vector to hold the dof indices on
    // the cell we visit at a time
    std::vector<types::global_dof_index> local_dof_indices;

    // prebuilt list to which component
    // a given dof on a cell
    // should go. note that we get into
    // trouble here if the shape
    // function is not primitive, since
    // then there is no single vector
    // component to which it
    // belongs. in this case, assign it
    // to the first vector component to
    // which it belongs
    std::vector<std::vector<unsigned int>> component_list(fe_collection.size());
    for (unsigned int f = 0; f < fe_collection.size(); ++f)
      {
        const FiniteElement<dim, spacedim> &fe = fe_collection[f];
        const unsigned int dofs_per_cell       = fe.n_dofs_per_cell();
        component_list[f].resize(dofs_per_cell);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          if (fe.is_primitive(i))
            component_list[f][i] =
              component_order[fe.system_to_component_index(i).first];
          else
            {
              const unsigned int comp =
                fe.get_nonzero_components(i).first_selected_component();

              // then associate this degree
              // of freedom with this
              // component
              component_list[f][i] = component_order[comp];
            }
      }

    // set up a map where for each
    // component the respective degrees
    // of freedom are collected.
    //
    // note that this map is sorted by
    // component but that within each
    // component it is NOT sorted by
    // dof index. note also that some
    // dof indices are entered
    // multiply, so we will have to
    // take care of that
    std::vector<std::vector<types::global_dof_index>> component_to_dof_map(
      fe_collection.n_components());
    for (CellIterator cell = start; cell != end; ++cell)
      {
        if (is_level_operation)
          {
            // we are dealing with mg dofs, skip foreign level cells:
            if (!cell->is_locally_owned_on_level())
              continue;
          }
        else
          {
            // we are dealing with active dofs, skip the loop if not locally
            // owned:
            if (!cell->is_locally_owned())
              continue;
          }

        if (is_level_operation)
          Assert(
            cell->level() == start->level(),
            ExcMessage(
              "Multigrid renumbering in compute_component_wise() needs to be applied to a single level!"));

        // on each cell: get dof indices
        // and insert them into the global
        // list using their component
        const types::fe_index fe_index = cell->active_fe_index();
        const unsigned int    dofs_per_cell =
          fe_collection[fe_index].n_dofs_per_cell();
        local_dof_indices.resize(dofs_per_cell);
        cell->get_active_or_mg_dof_indices(local_dof_indices);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          if (locally_owned_dofs.is_element(local_dof_indices[i]))
            component_to_dof_map[component_list[fe_index][i]].push_back(
              local_dof_indices[i]);
      }

    // now we've got all indices sorted
    // into buckets labeled by their
    // target component number. we've
    // only got to traverse this list
    // and assign the new indices
    //
    // however, we first want to sort
    // the indices entered into the
    // buckets to preserve the order
    // within each component and during
    // this also remove duplicate
    // entries
    //
    // note that we no longer have to
    // care about non-primitive shape
    // functions since the buckets
    // corresponding to the second and
    // following vector components of a
    // non-primitive FE will simply be
    // empty, everything being shoved
    // into the first one. The same
    // holds if several components were
    // joined into a single target.
    for (unsigned int component = 0; component < fe_collection.n_components();
         ++component)
      {
        std::sort(component_to_dof_map[component].begin(),
                  component_to_dof_map[component].end());
        component_to_dof_map[component].erase(
          std::unique(component_to_dof_map[component].begin(),
                      component_to_dof_map[component].end()),
          component_to_dof_map[component].end());
      }

    // calculate the number of locally owned
    // DoFs per bucket
    const unsigned int n_buckets = fe_collection.n_components();
    std::vector<types::global_dof_index> shifts(n_buckets);

    if (const parallel::TriangulationBase<dim, spacedim> *tria =
          (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &start->get_dof_handler().get_triangulation())))
      {
#ifdef DEAL_II_WITH_MPI
        std::vector<types::global_dof_index> local_dof_count(n_buckets);

        for (unsigned int c = 0; c < n_buckets; ++c)
          local_dof_count[c] = component_to_dof_map[c].size();

        std::vector<types::global_dof_index> prefix_dof_count(n_buckets);
        const int                            ierr = MPI_Exscan(
          local_dof_count.data(),
          prefix_dof_count.data(),
          n_buckets,
          Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>,
          MPI_SUM,
          tria->get_mpi_communicator());
        AssertThrowMPI(ierr);

        std::vector<types::global_dof_index> global_dof_count(n_buckets);
        Utilities::MPI::sum(local_dof_count,
                            tria->get_mpi_communicator(),
                            global_dof_count);

        // calculate shifts
        types::global_dof_index cumulated = 0;
        for (unsigned int c = 0; c < n_buckets; ++c)
          {
            shifts[c] = prefix_dof_count[c] + cumulated;
            cumulated += global_dof_count[c];
          }
#else
        (void)tria;
        DEAL_II_ASSERT_UNREACHABLE();
#endif
      }
    else
      {
        shifts[0] = 0;
        for (unsigned int c = 1; c < fe_collection.n_components(); ++c)
          shifts[c] = shifts[c - 1] + component_to_dof_map[c - 1].size();
      }



    // now concatenate all the
    // components in the order the user
    // desired to see
    types::global_dof_index next_free_index = 0;
    for (unsigned int component = 0; component < fe_collection.n_components();
         ++component)
      {
        next_free_index = shifts[component];

        for (const types::global_dof_index dof_index :
             component_to_dof_map[component])
          {
            Assert(locally_owned_dofs.index_within_set(dof_index) <
                     new_indices.size(),
                   ExcInternalError());
            new_indices[locally_owned_dofs.index_within_set(dof_index)] =
              next_free_index;

            ++next_free_index;
          }
      }

    return next_free_index;
  }



  template <int dim, int spacedim>
  void
  block_wise(DoFHandler<dim, spacedim> &dof_handler)
  {
    std::vector<types::global_dof_index> renumbering(
      dof_handler.n_locally_owned_dofs(), numbers::invalid_dof_index);

    const types::global_dof_index result = compute_block_wise<
      dim,
      spacedim,
      typename DoFHandler<dim, spacedim>::active_cell_iterator,
      typename DoFHandler<dim, spacedim>::level_cell_iterator>(
      renumbering, dof_handler.begin_active(), dof_handler.end(), false);
    if (result == 0)
      return;

    // verify that the last numbered
    // degree of freedom is either
    // equal to the number of degrees
    // of freedom in total (the
    // sequential case) or in the
    // distributed case at least
    // makes sense
    Assert((result == dof_handler.n_locally_owned_dofs()) ||
             ((dof_handler.n_locally_owned_dofs() < dof_handler.n_dofs()) &&
              (result <= dof_handler.n_dofs())),
           ExcInternalError());

    dof_handler.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  block_wise(DoFHandler<dim, spacedim> &dof_handler, const unsigned int level)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcDoFHandlerNotInitialized());

    std::vector<types::global_dof_index> renumbering(
      dof_handler.locally_owned_mg_dofs(level).n_elements(),
      numbers::invalid_dof_index);

    typename DoFHandler<dim, spacedim>::level_cell_iterator start =
      dof_handler.begin(level);
    typename DoFHandler<dim, spacedim>::level_cell_iterator end =
      dof_handler.end(level);

    const types::global_dof_index result = compute_block_wise<
      dim,
      spacedim,
      typename DoFHandler<dim, spacedim>::level_cell_iterator,
      typename DoFHandler<dim, spacedim>::level_cell_iterator>(renumbering,
                                                               start,
                                                               end,
                                                               true);
    (void)result;

    Assert(result == 0 || result == dof_handler.n_dofs(level),
           ExcInternalError());

    if (Utilities::MPI::max(renumbering.size(),
                            dof_handler.get_mpi_communicator()) > 0)
      dof_handler.renumber_dofs(level, renumbering);
  }



  template <int dim, int spacedim, class IteratorType, class EndIteratorType>
  types::global_dof_index
  compute_block_wise(std::vector<types::global_dof_index> &new_indices,
                     const IteratorType                   &start,
                     const EndIteratorType                &end,
                     const bool                            is_level_operation)
  {
    const hp::FECollection<dim, spacedim> &fe_collection =
      start->get_dof_handler().get_fe_collection();

    // do nothing if the FE has only
    // one component
    if (fe_collection.n_blocks() == 1)
      {
        new_indices.resize(0);
        return 0;
      }

    // Get a reference to the set of dofs. Note that we assume that all cells
    // are assumed to be on the same level, otherwise the operation doesn't make
    // much sense (we will assert this below).
    const IndexSet &locally_owned_dofs =
      is_level_operation ?
        start->get_dof_handler().locally_owned_mg_dofs(start->level()) :
        start->get_dof_handler().locally_owned_dofs();

    // vector to hold the dof indices on
    // the cell we visit at a time
    std::vector<types::global_dof_index> local_dof_indices;

    // prebuilt list to which block
    // a given dof on a cell
    // should go.
    std::vector<std::vector<types::global_dof_index>> block_list(
      fe_collection.size());
    for (unsigned int f = 0; f < fe_collection.size(); ++f)
      {
        const FiniteElement<dim, spacedim> &fe = fe_collection[f];
        block_list[f].resize(fe.n_dofs_per_cell());
        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          block_list[f][i] = fe.system_to_block_index(i).first;
      }

    // set up a map where for each
    // block the respective degrees
    // of freedom are collected.
    //
    // note that this map is sorted by
    // block but that within each
    // block it is NOT sorted by
    // dof index. note also that some
    // dof indices are entered
    // multiply, so we will have to
    // take care of that
    std::vector<std::vector<types::global_dof_index>> block_to_dof_map(
      fe_collection.n_blocks());
    for (IteratorType cell = start; cell != end; ++cell)
      {
        if (is_level_operation)
          {
            // we are dealing with mg dofs, skip foreign level cells:
            if (!cell->is_locally_owned_on_level())
              continue;
          }
        else
          {
            // we are dealing with active dofs, skip the loop if not locally
            // owned:
            if (!cell->is_locally_owned())
              continue;
          }

        if (is_level_operation)
          Assert(
            cell->level() == start->level(),
            ExcMessage(
              "Multigrid renumbering in compute_block_wise() needs to be applied to a single level!"));

        // on each cell: get dof indices
        // and insert them into the global
        // list using their component
        const types::fe_index fe_index = cell->active_fe_index();
        const unsigned int    dofs_per_cell =
          fe_collection[fe_index].n_dofs_per_cell();
        local_dof_indices.resize(dofs_per_cell);
        cell->get_active_or_mg_dof_indices(local_dof_indices);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          if (locally_owned_dofs.is_element(local_dof_indices[i]))
            block_to_dof_map[block_list[fe_index][i]].push_back(
              local_dof_indices[i]);
      }

    // now we've got all indices sorted
    // into buckets labeled by their
    // target block number. we've
    // only got to traverse this list
    // and assign the new indices
    //
    // however, we first want to sort
    // the indices entered into the
    // buckets to preserve the order
    // within each component and during
    // this also remove duplicate
    // entries
    for (unsigned int block = 0; block < fe_collection.n_blocks(); ++block)
      {
        std::sort(block_to_dof_map[block].begin(),
                  block_to_dof_map[block].end());
        block_to_dof_map[block].erase(
          std::unique(block_to_dof_map[block].begin(),
                      block_to_dof_map[block].end()),
          block_to_dof_map[block].end());
      }

    // calculate the number of locally owned
    // DoFs per bucket
    const unsigned int                   n_buckets = fe_collection.n_blocks();
    std::vector<types::global_dof_index> shifts(n_buckets);

    if (const parallel::TriangulationBase<dim, spacedim> *tria =
          (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &start->get_dof_handler().get_triangulation())))
      {
#ifdef DEAL_II_WITH_MPI
        std::vector<types::global_dof_index> local_dof_count(n_buckets);

        for (unsigned int c = 0; c < n_buckets; ++c)
          local_dof_count[c] = block_to_dof_map[c].size();

        std::vector<types::global_dof_index> prefix_dof_count(n_buckets);
        const int                            ierr = MPI_Exscan(
          local_dof_count.data(),
          prefix_dof_count.data(),
          n_buckets,
          Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>,
          MPI_SUM,
          tria->get_mpi_communicator());
        AssertThrowMPI(ierr);

        std::vector<types::global_dof_index> global_dof_count(n_buckets);
        Utilities::MPI::sum(local_dof_count,
                            tria->get_mpi_communicator(),
                            global_dof_count);

        // calculate shifts
        types::global_dof_index cumulated = 0;
        for (unsigned int c = 0; c < n_buckets; ++c)
          {
            shifts[c] = prefix_dof_count[c] + cumulated;
            cumulated += global_dof_count[c];
          }
#else
        (void)tria;
        DEAL_II_ASSERT_UNREACHABLE();
#endif
      }
    else
      {
        shifts[0] = 0;
        for (unsigned int c = 1; c < fe_collection.n_blocks(); ++c)
          shifts[c] = shifts[c - 1] + block_to_dof_map[c - 1].size();
      }



    // now concatenate all the
    // components in the order the user
    // desired to see
    types::global_dof_index next_free_index = 0;
    for (unsigned int block = 0; block < fe_collection.n_blocks(); ++block)
      {
        const typename std::vector<types::global_dof_index>::const_iterator
          begin_of_component = block_to_dof_map[block].begin(),
          end_of_component   = block_to_dof_map[block].end();

        next_free_index = shifts[block];

        for (typename std::vector<types::global_dof_index>::const_iterator
               dof_index = begin_of_component;
             dof_index != end_of_component;
             ++dof_index)
          {
            Assert(locally_owned_dofs.index_within_set(*dof_index) <
                     new_indices.size(),
                   ExcInternalError());
            new_indices[locally_owned_dofs.index_within_set(*dof_index)] =
              next_free_index++;
          }
      }

    return next_free_index;
  }



  namespace
  {
    // Helper function for DoFRenumbering::hierarchical(). This function
    // recurses into the given cell or, if that should be an active (terminal)
    // cell, renumbers DoF indices on it. The function starts renumbering with
    // 'next_free_dof_index' and returns the first still unused DoF index at the
    // end of its operation.
    template <int dim, typename CellIteratorType>
    types::global_dof_index
    compute_hierarchical_recursive(
      const types::global_dof_index         next_free_dof_offset,
      const types::global_dof_index         my_starting_index,
      const CellIteratorType               &cell,
      const IndexSet                       &locally_owned_dof_indices,
      std::vector<types::global_dof_index> &new_indices)
    {
      types::global_dof_index current_next_free_dof_offset =
        next_free_dof_offset;

      if (cell->has_children())
        {
          // recursion
          for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
               ++c)
            current_next_free_dof_offset =
              compute_hierarchical_recursive<dim>(current_next_free_dof_offset,
                                                  my_starting_index,
                                                  cell->child(c),
                                                  locally_owned_dof_indices,
                                                  new_indices);
        }
      else
        {
          // this is a terminal cell. we need to renumber its DoF indices. there
          // are now three cases to decide:
          // - this is a sequential triangulation: we can just go ahead and
          // number
          //   the DoFs in the order in which we encounter cells. in this case,
          //   all cells are actually locally owned
          // - if this is a parallel::distributed::Triangulation, then we only
          //   need to work on the locally owned cells since they contain
          //   all locally owned DoFs.
          // - if this is a parallel::shared::Triangulation, then the same
          // applies
          //
          // in all cases, each processor starts new indices so that we get
          // a consecutive numbering on each processor, and disjoint ownership
          // of the global range of DoF indices
          if (cell->is_locally_owned())
            {
              // first get the existing DoF indices
              const unsigned int dofs_per_cell =
                cell->get_fe().n_dofs_per_cell();
              std::vector<types::global_dof_index> local_dof_indices(
                dofs_per_cell);
              cell->get_dof_indices(local_dof_indices);

              // then loop over the existing DoF indices on this cell
              // and see whether it has already been re-numbered (it
              // may have been on a face or vertex to a neighboring
              // cell that we have encountered before). if not,
              // give it a new number and store that number in the
              // output array (new_indices)
              //
              // if this is a parallel triangulation and a DoF index is
              // not locally owned, then don't touch it. since
              // we don't actually *change* DoF indices (just record new
              // numbers in an array), we don't need to worry about
              // the decision whether a DoF is locally owned or not changing
              // as we progress in renumbering DoFs -- all adjacent cells
              // will always agree that a DoF is locally owned or not.
              // that said, the first cell to encounter a locally owned DoF
              // gets to number it, so the order in which we traverse cells
              // matters
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                if (locally_owned_dof_indices.is_element(local_dof_indices[i]))
                  {
                    // this is a locally owned DoF, assign new number if not
                    // assigned a number yet
                    const unsigned int idx =
                      locally_owned_dof_indices.index_within_set(
                        local_dof_indices[i]);
                    if (new_indices[idx] == numbers::invalid_dof_index)
                      {
                        new_indices[idx] =
                          my_starting_index + current_next_free_dof_offset;
                        ++current_next_free_dof_offset;
                      }
                  }
            }
        }

      return current_next_free_dof_offset;
    }
  } // namespace



  template <int dim, int spacedim>
  void
  hierarchical(DoFHandler<dim, spacedim> &dof_handler)
  {
    std::vector<types::global_dof_index> renumbering(
      dof_handler.n_locally_owned_dofs(), numbers::invalid_dof_index);

    types::global_dof_index next_free_dof_offset = 0;
    const IndexSet          locally_owned = dof_handler.locally_owned_dofs();

    // in the function we call recursively, we want to number DoFs so
    // that global cell zero starts with DoF zero, regardless of how
    // DoFs were previously numbered. to this end, we need to figure
    // out which DoF index the current processor should start with.
    //
    // if this is a sequential triangulation, then obviously the starting
    // index is zero. otherwise, make sure we get contiguous, successive
    // ranges on each processor. note that since the number of locally owned
    // DoFs is an invariant under renumbering, we can easily compute this
    // starting index by just accumulating over the number of locally owned
    // DoFs for all previous processes
    types::global_dof_index my_starting_index = 0;

    if (const parallel::TriangulationBase<dim, spacedim> *tria =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &dof_handler.get_triangulation()))
      {
#ifdef DEAL_II_WITH_MPI
        types::global_dof_index locally_owned_size =
          dof_handler.locally_owned_dofs().n_elements();
        const int ierr = MPI_Exscan(
          &locally_owned_size,
          &my_starting_index,
          1,
          Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>,
          MPI_SUM,
          tria->get_mpi_communicator());
        AssertThrowMPI(ierr);
#endif
      }

    if (const parallel::distributed::Triangulation<dim, spacedim> *tria =
          dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim>
                         *>(&dof_handler.get_triangulation()))
      {
#ifdef DEAL_II_WITH_P4EST
        // this is a distributed Triangulation. we need to traverse the coarse
        // cells in the order p4est does to match the z-order actually used
        // by p4est. this requires using the renumbering of coarse cells
        // we do before we hand things off to p4est
        for (unsigned int c = 0; c < tria->n_cells(0); ++c)
          {
            const unsigned int coarse_cell_index =
              tria->get_p4est_tree_to_coarse_cell_permutation()[c];

            const typename DoFHandler<dim, spacedim>::level_cell_iterator
              this_cell(tria, 0, coarse_cell_index, &dof_handler);

            next_free_dof_offset =
              compute_hierarchical_recursive<dim>(next_free_dof_offset,
                                                  my_starting_index,
                                                  this_cell,
                                                  locally_owned,
                                                  renumbering);
          }
#else
        DEAL_II_NOT_IMPLEMENTED();
#endif
      }
    else
      {
        // this is not a distributed Triangulation, so we can traverse coarse
        // cells in the normal order
        for (typename DoFHandler<dim, spacedim>::cell_iterator cell =
               dof_handler.begin(0);
             cell != dof_handler.end(0);
             ++cell)
          next_free_dof_offset =
            compute_hierarchical_recursive<dim>(next_free_dof_offset,
                                                my_starting_index,
                                                cell,
                                                locally_owned,
                                                renumbering);
      }

    // verify that the last numbered degree of freedom is either
    // equal to the number of degrees of freedom in total (the
    // sequential case) or in the distributed case at least
    // makes sense
    Assert((next_free_dof_offset == dof_handler.n_locally_owned_dofs()) ||
             ((dof_handler.n_locally_owned_dofs() < dof_handler.n_dofs()) &&
              (next_free_dof_offset <= dof_handler.n_dofs())),
           ExcInternalError());

    // make sure that all local DoFs got new numbers assigned
    Assert(std::find(renumbering.begin(),
                     renumbering.end(),
                     numbers::invalid_dof_index) == renumbering.end(),
           ExcInternalError());

    dof_handler.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  sort_selected_dofs_back(DoFHandler<dim, spacedim> &dof_handler,
                          const std::vector<bool>   &selected_dofs)
  {
    std::vector<types::global_dof_index> renumbering(
      dof_handler.n_dofs(), numbers::invalid_dof_index);
    compute_sort_selected_dofs_back(renumbering, dof_handler, selected_dofs);

    dof_handler.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  sort_selected_dofs_back(DoFHandler<dim, spacedim> &dof_handler,
                          const std::vector<bool>   &selected_dofs,
                          const unsigned int         level)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcDoFHandlerNotInitialized());

    std::vector<types::global_dof_index> renumbering(
      dof_handler.n_dofs(level), numbers::invalid_dof_index);
    compute_sort_selected_dofs_back(renumbering,
                                    dof_handler,
                                    selected_dofs,
                                    level);

    dof_handler.renumber_dofs(level, renumbering);
  }



  template <int dim, int spacedim>
  void
  compute_sort_selected_dofs_back(
    std::vector<types::global_dof_index> &new_indices,
    const DoFHandler<dim, spacedim>      &dof_handler,
    const std::vector<bool>              &selected_dofs)
  {
    const types::global_dof_index n_dofs = dof_handler.n_dofs();
    Assert(selected_dofs.size() == n_dofs,
           ExcDimensionMismatch(selected_dofs.size(), n_dofs));

    // re-sort the dofs according to
    // their selection state
    Assert(new_indices.size() == n_dofs,
           ExcDimensionMismatch(new_indices.size(), n_dofs));

    const types::global_dof_index n_selected_dofs =
      std::count(selected_dofs.begin(), selected_dofs.end(), false);

    types::global_dof_index next_unselected = 0;
    types::global_dof_index next_selected   = n_selected_dofs;
    for (types::global_dof_index i = 0; i < n_dofs; ++i)
      if (selected_dofs[i] == false)
        {
          new_indices[i] = next_unselected;
          ++next_unselected;
        }
      else
        {
          new_indices[i] = next_selected;
          ++next_selected;
        }
    Assert(next_unselected == n_selected_dofs, ExcInternalError());
    Assert(next_selected == n_dofs, ExcInternalError());
  }



  template <int dim, int spacedim>
  void
  compute_sort_selected_dofs_back(
    std::vector<types::global_dof_index> &new_indices,
    const DoFHandler<dim, spacedim>      &dof_handler,
    const std::vector<bool>              &selected_dofs,
    const unsigned int                    level)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcDoFHandlerNotInitialized());

    const unsigned int n_dofs = dof_handler.n_dofs(level);
    Assert(selected_dofs.size() == n_dofs,
           ExcDimensionMismatch(selected_dofs.size(), n_dofs));

    // re-sort the dofs according to
    // their selection state
    Assert(new_indices.size() == n_dofs,
           ExcDimensionMismatch(new_indices.size(), n_dofs));

    const unsigned int n_selected_dofs =
      std::count(selected_dofs.begin(), selected_dofs.end(), false);

    unsigned int next_unselected = 0;
    unsigned int next_selected   = n_selected_dofs;
    for (unsigned int i = 0; i < n_dofs; ++i)
      if (selected_dofs[i] == false)
        {
          new_indices[i] = next_unselected;
          ++next_unselected;
        }
      else
        {
          new_indices[i] = next_selected;
          ++next_selected;
        }
    Assert(next_unselected == n_selected_dofs, ExcInternalError());
    Assert(next_selected == n_dofs, ExcInternalError());
  }



  template <int dim, int spacedim>
  void
  cell_wise(
    DoFHandler<dim, spacedim> &dof,
    const std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
      &cells)
  {
    std::vector<types::global_dof_index> renumbering(
      dof.n_locally_owned_dofs());
    std::vector<types::global_dof_index> reverse(dof.n_locally_owned_dofs());
    compute_cell_wise(renumbering, reverse, dof, cells);

    dof.renumber_dofs(renumbering);
  }


  template <int dim, int spacedim>
  void
  compute_cell_wise(
    std::vector<types::global_dof_index> &new_indices,
    std::vector<types::global_dof_index> &reverse,
    const DoFHandler<dim, spacedim>      &dof,
    const typename std::vector<
      typename DoFHandler<dim, spacedim>::active_cell_iterator> &cells)
  {
    if (const parallel::TriangulationBase<dim, spacedim> *p =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &dof.get_triangulation()))
      {
        AssertDimension(cells.size(), p->n_locally_owned_active_cells());
      }
    else
      {
        AssertDimension(cells.size(), dof.get_triangulation().n_active_cells());
      }

    const auto n_owned_dofs = dof.n_locally_owned_dofs();

    // Actually, we compute the inverse of the reordering vector, called reverse
    // here. Later, its inverse is computed into new_indices, which is the
    // return argument.

    AssertDimension(new_indices.size(), n_owned_dofs);
    AssertDimension(reverse.size(), n_owned_dofs);

    // For continuous elements, we must make sure, that each dof is reordered
    // only once.
    std::vector<bool>                    already_sorted(n_owned_dofs, false);
    std::vector<types::global_dof_index> cell_dofs;

    const auto &owned_dofs = dof.locally_owned_dofs();

    unsigned int index = 0;

    for (const auto &cell : cells)
      {
        // Determine the number of dofs on this cell and reinit the
        // vector storing these numbers.
        const unsigned int n_cell_dofs = cell->get_fe().n_dofs_per_cell();
        cell_dofs.resize(n_cell_dofs);

        cell->get_active_or_mg_dof_indices(cell_dofs);

        // Sort here to make sure that degrees of freedom inside a single cell
        // are in the same order after renumbering.
        std::sort(cell_dofs.begin(), cell_dofs.end());

        for (const auto dof : cell_dofs)
          {
            const auto local_dof = owned_dofs.index_within_set(dof);
            if (local_dof != numbers::invalid_dof_index &&
                !already_sorted[local_dof])
              {
                already_sorted[local_dof] = true;
                reverse[index++]          = local_dof;
              }
          }
      }
    Assert(index == n_owned_dofs,
           ExcMessage(
             "Traversing over the given set of cells did not cover all "
             "degrees of freedom in the DoFHandler. Does the set of cells "
             "not include all active cells?"));

    for (types::global_dof_index i = 0; i < reverse.size(); ++i)
      new_indices[reverse[i]] = owned_dofs.nth_index_in_set(i);
  }



  template <int dim, int spacedim>
  void
  cell_wise(DoFHandler<dim, spacedim> &dof,
            const unsigned int         level,
            const typename std::vector<
              typename DoFHandler<dim, spacedim>::level_cell_iterator> &cells)
  {
    Assert(dof.n_dofs(level) != numbers::invalid_dof_index,
           ExcDoFHandlerNotInitialized());

    std::vector<types::global_dof_index> renumbering(dof.n_dofs(level));
    std::vector<types::global_dof_index> reverse(dof.n_dofs(level));

    compute_cell_wise(renumbering, reverse, dof, level, cells);
    dof.renumber_dofs(level, renumbering);
  }



  template <int dim, int spacedim>
  void
  compute_cell_wise(
    std::vector<types::global_dof_index> &new_order,
    std::vector<types::global_dof_index> &reverse,
    const DoFHandler<dim, spacedim>      &dof,
    const unsigned int                    level,
    const typename std::vector<
      typename DoFHandler<dim, spacedim>::level_cell_iterator> &cells)
  {
    Assert(cells.size() == dof.get_triangulation().n_cells(level),
           ExcDimensionMismatch(cells.size(),
                                dof.get_triangulation().n_cells(level)));
    Assert(new_order.size() == dof.n_dofs(level),
           ExcDimensionMismatch(new_order.size(), dof.n_dofs(level)));
    Assert(reverse.size() == dof.n_dofs(level),
           ExcDimensionMismatch(reverse.size(), dof.n_dofs(level)));

    unsigned int n_global_dofs = dof.n_dofs(level);
    unsigned int n_cell_dofs   = dof.get_fe().n_dofs_per_cell();

    std::vector<bool>                    already_sorted(n_global_dofs, false);
    std::vector<types::global_dof_index> cell_dofs(n_cell_dofs);

    unsigned int global_index = 0;

    for (const auto &cell : cells)
      {
        Assert(cell->level() == static_cast<int>(level), ExcInternalError());

        cell->get_active_or_mg_dof_indices(cell_dofs);
        std::sort(cell_dofs.begin(), cell_dofs.end());

        for (unsigned int i = 0; i < n_cell_dofs; ++i)
          {
            if (!already_sorted[cell_dofs[i]])
              {
                already_sorted[cell_dofs[i]] = true;
                reverse[global_index++]      = cell_dofs[i];
              }
          }
      }
    Assert(global_index == n_global_dofs,
           ExcMessage(
             "Traversing over the given set of cells did not cover all "
             "degrees of freedom in the DoFHandler. Does the set of cells "
             "not include all cells of the specified level?"));

    for (types::global_dof_index i = 0; i < new_order.size(); ++i)
      new_order[reverse[i]] = i;
  }



  template <int dim, int spacedim>
  void
  downstream(DoFHandler<dim, spacedim> &dof,
             const Tensor<1, spacedim> &direction,
             const bool                 dof_wise_renumbering)
  {
    std::vector<types::global_dof_index> renumbering(dof.n_dofs());
    std::vector<types::global_dof_index> reverse(dof.n_dofs());
    compute_downstream(
      renumbering, reverse, dof, direction, dof_wise_renumbering);

    dof.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  compute_downstream(std::vector<types::global_dof_index> &new_indices,
                     std::vector<types::global_dof_index> &reverse,
                     const DoFHandler<dim, spacedim>      &dof,
                     const Tensor<1, spacedim>            &direction,
                     const bool                            dof_wise_renumbering)
  {
    Assert((dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
              &dof.get_triangulation()) == nullptr),
           ExcNotImplemented());

    if (dof_wise_renumbering == false)
      {
        std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
          ordered_cells;
        ordered_cells.reserve(dof.get_triangulation().n_active_cells());
        const CompareDownstream<
          typename DoFHandler<dim, spacedim>::active_cell_iterator,
          spacedim>
          comparator(direction);

        for (const auto &cell : dof.active_cell_iterators())
          ordered_cells.push_back(cell);

        std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);

        compute_cell_wise(new_indices, reverse, dof, ordered_cells);
      }
    else
      {
        // similar code as for
        // DoFTools::map_dofs_to_support_points, but
        // need to do this for general DoFHandler<dim, spacedim> classes and
        // want to be able to sort the result
        // (otherwise, could use something like
        // DoFTools::map_support_points_to_dofs)
        const unsigned int n_dofs = dof.n_dofs();
        std::vector<std::pair<Point<spacedim>, unsigned int>>
          support_point_list(n_dofs);

        const hp::FECollection<dim> &fe_collection = dof.get_fe_collection();
        Assert(fe_collection[0].has_support_points(),
               typename FiniteElement<dim>::ExcFEHasNoSupportPoints());
        hp::QCollection<dim> quadrature_collection;
        for (unsigned int comp = 0; comp < fe_collection.size(); ++comp)
          {
            Assert(fe_collection[comp].has_support_points(),
                   typename FiniteElement<dim>::ExcFEHasNoSupportPoints());
            quadrature_collection.push_back(
              Quadrature<dim>(fe_collection[comp].get_unit_support_points()));
          }
        hp::FEValues<dim, spacedim> hp_fe_values(fe_collection,
                                                 quadrature_collection,
                                                 update_quadrature_points);

        std::vector<bool> already_touched(n_dofs, false);

        std::vector<types::global_dof_index> local_dof_indices;

        for (const auto &cell : dof.active_cell_iterators())
          {
            const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
            local_dof_indices.resize(dofs_per_cell);
            hp_fe_values.reinit(cell);
            const FEValues<dim> &fe_values =
              hp_fe_values.get_present_fe_values();
            cell->get_active_or_mg_dof_indices(local_dof_indices);
            const std::vector<Point<spacedim>> &points =
              fe_values.get_quadrature_points();
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              if (!already_touched[local_dof_indices[i]])
                {
                  support_point_list[local_dof_indices[i]].first = points[i];
                  support_point_list[local_dof_indices[i]].second =
                    local_dof_indices[i];
                  already_touched[local_dof_indices[i]] = true;
                }
          }

        ComparePointwiseDownstream<spacedim> comparator(direction);
        std::sort(support_point_list.begin(),
                  support_point_list.end(),
                  comparator);
        for (types::global_dof_index i = 0; i < n_dofs; ++i)
          new_indices[support_point_list[i].second] = i;
      }
  }



  template <int dim, int spacedim>
  void
  downstream(DoFHandler<dim, spacedim> &dof,
             const unsigned int         level,
             const Tensor<1, spacedim> &direction,
             const bool                 dof_wise_renumbering)
  {
    std::vector<types::global_dof_index> renumbering(dof.n_dofs(level));
    std::vector<types::global_dof_index> reverse(dof.n_dofs(level));
    compute_downstream(
      renumbering, reverse, dof, level, direction, dof_wise_renumbering);

    dof.renumber_dofs(level, renumbering);
  }



  template <int dim, int spacedim>
  void
  compute_downstream(std::vector<types::global_dof_index> &new_indices,
                     std::vector<types::global_dof_index> &reverse,
                     const DoFHandler<dim, spacedim>      &dof,
                     const unsigned int                    level,
                     const Tensor<1, spacedim>            &direction,
                     const bool                            dof_wise_renumbering)
  {
    if (dof_wise_renumbering == false)
      {
        std::vector<typename DoFHandler<dim, spacedim>::level_cell_iterator>
          ordered_cells;
        ordered_cells.reserve(dof.get_triangulation().n_cells(level));
        const CompareDownstream<
          typename DoFHandler<dim, spacedim>::level_cell_iterator,
          spacedim>
          comparator(direction);

        typename DoFHandler<dim, spacedim>::level_cell_iterator p =
          dof.begin(level);
        typename DoFHandler<dim, spacedim>::level_cell_iterator end =
          dof.end(level);

        while (p != end)
          {
            ordered_cells.push_back(p);
            ++p;
          }
        std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);

        compute_cell_wise(new_indices, reverse, dof, level, ordered_cells);
      }
    else
      {
        Assert(dof.get_fe().has_support_points(),
               typename FiniteElement<dim>::ExcFEHasNoSupportPoints());
        const unsigned int n_dofs = dof.n_dofs(level);
        std::vector<std::pair<Point<spacedim>, unsigned int>>
          support_point_list(n_dofs);

        const Quadrature<dim>   q_dummy(dof.get_fe().get_unit_support_points());
        FEValues<dim, spacedim> fe_values(dof.get_fe(),
                                          q_dummy,
                                          update_quadrature_points);

        std::vector<bool> already_touched(dof.n_dofs(), false);

        const unsigned int dofs_per_cell = dof.get_fe().n_dofs_per_cell();
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        typename DoFHandler<dim, spacedim>::level_cell_iterator begin =
          dof.begin(level);
        typename DoFHandler<dim, spacedim>::level_cell_iterator end =
          dof.end(level);
        for (; begin != end; ++begin)
          {
            const typename Triangulation<dim, spacedim>::cell_iterator
              &begin_tria = begin;
            begin->get_active_or_mg_dof_indices(local_dof_indices);
            fe_values.reinit(begin_tria);
            const std::vector<Point<spacedim>> &points =
              fe_values.get_quadrature_points();
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              if (!already_touched[local_dof_indices[i]])
                {
                  support_point_list[local_dof_indices[i]].first = points[i];
                  support_point_list[local_dof_indices[i]].second =
                    local_dof_indices[i];
                  already_touched[local_dof_indices[i]] = true;
                }
          }

        ComparePointwiseDownstream<spacedim> comparator(direction);
        std::sort(support_point_list.begin(),
                  support_point_list.end(),
                  comparator);
        for (types::global_dof_index i = 0; i < n_dofs; ++i)
          new_indices[support_point_list[i].second] = i;
      }
  }



  /**
   * Provide comparator for DoFCellAccessors
   */
  namespace internal
  {
    template <int dim>
    struct ClockCells
    {
      /**
       * Center of rotation.
       */
      const Point<dim> &center;
      /**
       * Revert sorting order.
       */
      bool counter;

      /**
       * Constructor.
       */
      ClockCells(const Point<dim> &center, bool counter)
        : center(center)
        , counter(counter)
      {}

      /**
       * Comparison operator
       */
      template <class DHCellIterator>
      bool
      operator()(const DHCellIterator &c1, const DHCellIterator &c2) const
      {
        // dispatch to
        // dimension-dependent functions
        return compare(c1, c2, std::integral_constant<int, dim>());
      }

    private:
      /**
       * Comparison operator for dim>=2
       */
      template <class DHCellIterator, int xdim>
      bool
      compare(const DHCellIterator &c1,
              const DHCellIterator &c2,
              std::integral_constant<int, xdim>) const
      {
        const Tensor<1, dim> v1 = c1->center() - center;
        const Tensor<1, dim> v2 = c2->center() - center;
        const double         s1 = std::atan2(v1[0], v1[1]);
        const double         s2 = std::atan2(v2[0], v2[1]);
        return (counter ? (s1 > s2) : (s2 > s1));
      }


      /**
       * Comparison operator for dim==1
       * where this function makes no sense
       */
      template <class DHCellIterator>
      bool
      compare(const DHCellIterator &,
              const DHCellIterator &,
              std::integral_constant<int, 1>) const
      {
        Assert(dim >= 2,
               ExcMessage("This operation only makes sense for dim>=2."));
        return false;
      }
    };
  } // namespace internal



  template <int dim, int spacedim>
  void
  clockwise_dg(DoFHandler<dim, spacedim> &dof,
               const Point<spacedim>     &center,
               const bool                 counter)
  {
    std::vector<types::global_dof_index> renumbering(dof.n_dofs());
    compute_clockwise_dg(renumbering, dof, center, counter);

    dof.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  compute_clockwise_dg(std::vector<types::global_dof_index> &new_indices,
                       const DoFHandler<dim, spacedim>      &dof,
                       const Point<spacedim>                &center,
                       const bool                            counter)
  {
    std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
      ordered_cells;
    ordered_cells.reserve(dof.get_triangulation().n_active_cells());
    internal::ClockCells<spacedim> comparator(center, counter);

    for (const auto &cell : dof.active_cell_iterators())
      ordered_cells.push_back(cell);

    std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);

    std::vector<types::global_dof_index> reverse(new_indices.size());
    compute_cell_wise(new_indices, reverse, dof, ordered_cells);
  }



  template <int dim, int spacedim>
  void
  clockwise_dg(DoFHandler<dim, spacedim> &dof,
               const unsigned int         level,
               const Point<spacedim>     &center,
               const bool                 counter)
  {
    std::vector<typename DoFHandler<dim, spacedim>::level_cell_iterator>
      ordered_cells;
    ordered_cells.reserve(dof.get_triangulation().n_active_cells());
    internal::ClockCells<spacedim> comparator(center, counter);

    typename DoFHandler<dim, spacedim>::level_cell_iterator p =
      dof.begin(level);
    typename DoFHandler<dim, spacedim>::level_cell_iterator end =
      dof.end(level);

    while (p != end)
      {
        ordered_cells.push_back(p);
        ++p;
      }
    std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);

    cell_wise(dof, level, ordered_cells);
  }



  template <int dim, int spacedim>
  void
  random(DoFHandler<dim, spacedim> &dof_handler)
  {
    std::vector<types::global_dof_index> renumbering(
      dof_handler.n_dofs(), numbers::invalid_dof_index);
    compute_random(renumbering, dof_handler);

    dof_handler.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  random(DoFHandler<dim, spacedim> &dof_handler, const unsigned int level)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcDoFHandlerNotInitialized());

    std::vector<types::global_dof_index> renumbering(
      dof_handler.locally_owned_mg_dofs(level).n_elements(),
      numbers::invalid_dof_index);

    compute_random(renumbering, dof_handler, level);

    dof_handler.renumber_dofs(level, renumbering);
  }



  template <int dim, int spacedim>
  void
  compute_random(std::vector<types::global_dof_index> &new_indices,
                 const DoFHandler<dim, spacedim>      &dof_handler)
  {
    const types::global_dof_index n_dofs = dof_handler.n_dofs();
    Assert(new_indices.size() == n_dofs,
           ExcDimensionMismatch(new_indices.size(), n_dofs));

    std::iota(new_indices.begin(),
              new_indices.end(),
              types::global_dof_index(0));

    // shuffle the elements; the following is essentially std::shuffle (which
    // is new in C++11) but with a boost URNG
    // we could use std::mt19937 here but doing so results in compiler-dependent
    // output
    ::boost::mt19937 random_number_generator;
    for (unsigned int i = 1; i < n_dofs; ++i)
      {
        // get a random number between 0 and i (inclusive)
        const unsigned int j =
          ::boost::random::uniform_int_distribution<>(0, i)(
            random_number_generator);

        // if possible, swap the elements
        if (i != j)
          std::swap(new_indices[i], new_indices[j]);
      }
  }



  template <int dim, int spacedim>
  void
  compute_random(std::vector<types::global_dof_index> &new_indices,
                 const DoFHandler<dim, spacedim>      &dof_handler,
                 const unsigned int                    level)
  {
    const types::global_dof_index n_dofs = dof_handler.n_dofs(level);
    Assert(new_indices.size() == n_dofs,
           ExcDimensionMismatch(new_indices.size(), n_dofs));

    std::iota(new_indices.begin(),
              new_indices.end(),
              types::global_dof_index(0));

    // shuffle the elements; the following is essentially std::shuffle (which
    // is new in C++11) but with a boost URNG
    // we could use std::mt19937 here but doing so results in
    // compiler-dependent output
    ::boost::mt19937 random_number_generator;
    for (unsigned int i = 1; i < n_dofs; ++i)
      {
        // get a random number between 0 and i (inclusive)
        const unsigned int j =
          ::boost::random::uniform_int_distribution<>(0, i)(
            random_number_generator);

        // if possible, swap the elements
        if (i != j)
          std::swap(new_indices[i], new_indices[j]);
      }
  }



  template <int dim, int spacedim>
  void
  subdomain_wise(DoFHandler<dim, spacedim> &dof_handler)
  {
    Assert(
      (!dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &dof_handler.get_triangulation())),
      ExcMessage(
        "Parallel triangulations are already enumerated according to their MPI process id."));

    std::vector<types::global_dof_index> renumbering(
      dof_handler.n_dofs(), numbers::invalid_dof_index);
    compute_subdomain_wise(renumbering, dof_handler);

    dof_handler.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  compute_subdomain_wise(std::vector<types::global_dof_index> &new_dof_indices,
                         const DoFHandler<dim, spacedim>      &dof_handler)
  {
    const types::global_dof_index n_dofs = dof_handler.n_dofs();
    Assert(new_dof_indices.size() == n_dofs,
           ExcDimensionMismatch(new_dof_indices.size(), n_dofs));

    // first get the association of each dof
    // with a subdomain and determine the total
    // number of subdomain ids used
    std::vector<types::subdomain_id> subdomain_association(n_dofs);
    DoFTools::get_subdomain_association(dof_handler, subdomain_association);
    const unsigned int n_subdomains =
      *std::max_element(subdomain_association.begin(),
                        subdomain_association.end()) +
      1;

    // then renumber the subdomains by first
    // looking at those belonging to subdomain
    // 0, then those of subdomain 1, etc. note
    // that the algorithm is stable, i.e. if
    // two dofs i,j have i<j and belong to the
    // same subdomain, then they will be in
    // this order also after reordering
    std::fill(new_dof_indices.begin(),
              new_dof_indices.end(),
              numbers::invalid_dof_index);
    types::global_dof_index next_free_index = 0;
    for (types::subdomain_id subdomain = 0; subdomain < n_subdomains;
         ++subdomain)
      for (types::global_dof_index i = 0; i < n_dofs; ++i)
        if (subdomain_association[i] == subdomain)
          {
            Assert(new_dof_indices[i] == numbers::invalid_dof_index,
                   ExcInternalError());
            new_dof_indices[i] = next_free_index;
            ++next_free_index;
          }

    // we should have numbered all dofs
    Assert(next_free_index == n_dofs, ExcInternalError());
    Assert(std::find(new_dof_indices.begin(),
                     new_dof_indices.end(),
                     numbers::invalid_dof_index) == new_dof_indices.end(),
           ExcInternalError());
  }



  template <int dim, int spacedim>
  void
  support_point_wise(DoFHandler<dim, spacedim> &dof_handler)
  {
    std::vector<types::global_dof_index> renumbering(
      dof_handler.n_locally_owned_dofs(), numbers::invalid_dof_index);
    compute_support_point_wise(renumbering, dof_handler);

    // If there is only one component then there is nothing to do, so check
    // first:
    if (Utilities::MPI::max(renumbering.size(),
                            dof_handler.get_mpi_communicator()) > 0)
      dof_handler.renumber_dofs(renumbering);
  }



  template <int dim, int spacedim>
  void
  compute_support_point_wise(
    std::vector<types::global_dof_index> &new_dof_indices,
    const DoFHandler<dim, spacedim>      &dof_handler)
  {
    const types::global_dof_index n_dofs = dof_handler.n_locally_owned_dofs();
    Assert(new_dof_indices.size() == n_dofs,
           ExcDimensionMismatch(new_dof_indices.size(), n_dofs));

    // This renumbering occurs in three steps:
    // 1. Compute the component-wise renumbering so that all DoFs of component
    //    i are less than the DoFs of component i + 1.
    // 2. Compute a second renumbering component_to_nodal in which the
    //    renumbering is now, for two components, [u0, v0, u1, v1, ...]: i.e.,
    //    DoFs are first sorted by component and then by support point.
    // 3. Compose the two renumberings to obtain the final result.

    // Step 1:
    std::vector<types::global_dof_index> component_renumbering(
      n_dofs, numbers::invalid_dof_index);
    compute_component_wise<dim, spacedim>(component_renumbering,
                                          dof_handler.begin_active(),
                                          dof_handler.end(),
                                          std::vector<unsigned int>(),
                                          false);

    if constexpr (running_in_debug_mode())
      {
        {
          const std::vector<types::global_dof_index> dofs_per_component =
            DoFTools::count_dofs_per_fe_component(dof_handler, true);
          for (const auto &dpc : dofs_per_component)
            Assert(dofs_per_component[0] == dpc, ExcNotImplemented());
        }
      }
    const unsigned int n_components =
      dof_handler.get_fe_collection().n_components();
    Assert(dof_handler.n_dofs() % n_components == 0, ExcInternalError());
    const types::global_dof_index dofs_per_component =
      dof_handler.n_dofs() / n_components;
    const types::global_dof_index local_dofs_per_component =
      dof_handler.n_locally_owned_dofs() / n_components;

    // At this point we have no more communication to do - simplify things by
    // returning early if possible
    if (component_renumbering.empty())
      {
        new_dof_indices.resize(0);
        return;
      }
    std::fill(new_dof_indices.begin(),
              new_dof_indices.end(),
              numbers::invalid_dof_index);
    // This index set equals what dof_handler.locally_owned_dofs() would be if
    // we executed the componentwise renumbering.
    IndexSet component_renumbered_dofs(dof_handler.n_dofs());
    // DoFs in each component are now consecutive, which IndexSet::add_indices()
    // can exploit by avoiding calls to sort. Make use of that by adding DoFs
    // one component at a time:
    std::vector<types::global_dof_index> component_dofs(
      local_dofs_per_component);
    for (unsigned int component = 0; component < n_components; ++component)
      {
        for (std::size_t i = 0; i < local_dofs_per_component; ++i)
          component_dofs[i] =
            component_renumbering[n_components * i + component];
        component_renumbered_dofs.add_indices(component_dofs.begin(),
                                              component_dofs.end());
      }
    component_renumbered_dofs.compress();
    if constexpr (running_in_debug_mode())
      {
        {
          IndexSet component_renumbered_dofs2(dof_handler.n_dofs());
          component_renumbered_dofs2.add_indices(component_renumbering.begin(),
                                                 component_renumbering.end());
          Assert(component_renumbered_dofs2 == component_renumbered_dofs,
                 ExcInternalError());
        }
      }
    for (const FiniteElement<dim, spacedim> &fe :
         dof_handler.get_fe_collection())
      {
        AssertThrow(fe.dofs_per_cell == 0 || fe.has_support_points(),
                    ExcNotImplemented());
        for (unsigned int i = 0; i < fe.n_base_elements(); ++i)
          AssertThrow(
            fe.base_element(0).get_unit_support_points() ==
              fe.base_element(i).get_unit_support_points(),
            ExcMessage(
              "All base elements should have the same support points."));
      }

    std::vector<types::global_dof_index> component_to_nodal(
      dof_handler.n_locally_owned_dofs(), numbers::invalid_dof_index);

    // Step 2:
    std::vector<types::global_dof_index> cell_dofs;
    std::vector<types::global_dof_index> component_renumbered_cell_dofs;
    const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
    // Reuse the original index space for the new renumbering: it is the right
    // size and is contiguous on the current processor
    auto next_dof_it = locally_owned_dofs.begin();
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          const FiniteElement<dim, spacedim> &fe = cell->get_fe();
          cell_dofs.resize(fe.dofs_per_cell);
          component_renumbered_cell_dofs.resize(fe.dofs_per_cell);
          cell->get_dof_indices(cell_dofs);
          // Apply the component renumbering while skipping any ghost dofs. This
          // algorithm assumes that all locally owned DoFs before the component
          // renumbering are still locally owned afterwards (just with a new
          // index).
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            {
              if (locally_owned_dofs.is_element(cell_dofs[i]))
                {
                  const auto local_index =
                    locally_owned_dofs.index_within_set(cell_dofs[i]);
                  component_renumbered_cell_dofs[i] =
                    component_renumbering[local_index];
                }
              else
                {
                  component_renumbered_cell_dofs[i] =
                    numbers::invalid_dof_index;
                }
            }

          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            {
              if (fe.system_to_component_index(i).first == 0 &&
                  component_renumbered_dofs.is_element(
                    component_renumbered_cell_dofs[i]))
                {
                  for (unsigned int component = 0;
                       component < fe.n_components();
                       ++component)
                    {
                      // Since we are indexing in an odd way here it is much
                      // simpler to compute the permutation separately and
                      // combine it at the end instead of doing both at once
                      const auto local_index =
                        component_renumbered_dofs.index_within_set(
                          component_renumbered_cell_dofs[i] +
                          dofs_per_component * component);

                      if (component_to_nodal[local_index] ==
                          numbers::invalid_dof_index)
                        {
                          component_to_nodal[local_index] = *next_dof_it;
                          ++next_dof_it;
                        }
                    }
                }
            }
        }

    // Step 3:
    for (std::size_t i = 0; i < dof_handler.n_locally_owned_dofs(); ++i)
      {
        const auto local_index =
          component_renumbered_dofs.index_within_set(component_renumbering[i]);
        new_dof_indices[i] = component_to_nodal[local_index];
      }
  }



  template <int dim>
  void
  lexicographic(DoFHandler<dim> &dof_handler, const double tolerance)
  {
    std::vector<types::global_dof_index> renumbering;
    compute_lexicographic(renumbering, dof_handler, tolerance);
    dof_handler.renumber_dofs(renumbering);
  }


  template <int dim>
  void
  compute_lexicographic(std::vector<types::global_dof_index> &new_dof_indices,
                        const DoFHandler<dim>                &dof_handler,
                        const double                          tolerance)
  {
    Assert(dynamic_cast<const parallel::distributed::Triangulation<dim> *>(
             &dof_handler.get_triangulation()) == nullptr,
           ExcMessage(
             "Lexicographic renumbering is not implemented for distributed "
             "triangulations."));

    std::map<types::global_dof_index, Point<dim>> dof_location_map;
    DoFTools::map_dofs_to_support_points(MappingQ1<dim>(),
                                         dof_handler,
                                         dof_location_map);
    std::vector<std::pair<types::global_dof_index, Point<dim>>>
      dof_location_vector;

    dof_location_vector.reserve(dof_location_map.size());
    for (const auto &s : dof_location_map)
      dof_location_vector.push_back(s);

    std::sort(
      dof_location_vector.begin(),
      dof_location_vector.end(),
      [&](const std::pair<types::global_dof_index, Point<dim>> &p1,
          const std::pair<types::global_dof_index, Point<dim>> &p2) -> bool {
        for (int i = dim - 1; i >= 0; --i)
          {
            const double diff = p1.second(i) - p2.second(i);
            // Check if p1 is significantly smaller than p2 in this dimension
            if (diff < -tolerance)
              return true;
            // Check if p1 is significantly larger than p2 in this dimension
            if (diff > tolerance)
              return false;
            // Otherwise, coordinates are considered equal within tolerance for
            // this dimension, continue to the next dimension.
          }
        // All coordinates are within tolerance, points are considered
        // equivalent. Use the original DoF index as a tie-breaker to preserve
        // relative order.
        return p1.first < p2.first;
      });

    new_dof_indices.resize(dof_location_vector.size(),
                           numbers::invalid_dof_index);

    for (types::global_dof_index dof = 0; dof < dof_location_vector.size();
         ++dof)
      new_dof_indices[dof_location_vector[dof].first] = dof;
  }



  template <int dim,
            int spacedim,
            typename Number,
            typename VectorizedArrayType>
  void
  matrix_free_data_locality(
    DoFHandler<dim, spacedim>                          &dof_handler,
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free)
  {
    const std::vector<types::global_dof_index> new_global_numbers =
      compute_matrix_free_data_locality(dof_handler, matrix_free);
    if (matrix_free.get_mg_level() == dealii::numbers::invalid_unsigned_int)
      dof_handler.renumber_dofs(new_global_numbers);
    else
      dof_handler.renumber_dofs(matrix_free.get_mg_level(), new_global_numbers);
  }



  template <int dim, int spacedim, typename Number, typename AdditionalDataType>
  void
  matrix_free_data_locality(DoFHandler<dim, spacedim>       &dof_handler,
                            const AffineConstraints<Number> &constraints,
                            const AdditionalDataType        &matrix_free_data)
  {
    const std::vector<types::global_dof_index> new_global_numbers =
      compute_matrix_free_data_locality(dof_handler,
                                        constraints,
                                        matrix_free_data);
    if (matrix_free_data.mg_level == dealii::numbers::invalid_unsigned_int)
      dof_handler.renumber_dofs(new_global_numbers);
    else
      dof_handler.renumber_dofs(matrix_free_data.mg_level, new_global_numbers);
  }



  template <int dim, int spacedim, typename Number, typename AdditionalDataType>
  std::vector<types::global_dof_index>
  compute_matrix_free_data_locality(
    const DoFHandler<dim, spacedim> &dof_handler,
    const AffineConstraints<Number> &constraints,
    const AdditionalDataType        &matrix_free_data)
  {
    AdditionalDataType my_mf_data    = matrix_free_data;
    my_mf_data.initialize_mapping    = false;
    my_mf_data.tasks_parallel_scheme = AdditionalDataType::none;

    typename AdditionalDataType::MatrixFreeType separate_matrix_free;
    separate_matrix_free.reinit(dealii::MappingQ1<dim>(),
                                dof_handler,
                                constraints,
                                dealii::QGauss<1>(2),
                                my_mf_data);
    return compute_matrix_free_data_locality(dof_handler, separate_matrix_free);
  }



  // Implementation details for matrix-free renumbering
  namespace
  {
    // Compute a vector of lists with number of unknowns of the same category
    // in terms of influence from other MPI processes, starting with unknowns
    // touched only by the local process and finally a new set of indices
    // going to different MPI neighbors. Later passes of the algorithm will
    // re-order unknowns within each of these sets.
    std::vector<std::vector<unsigned int>>
    group_dofs_by_rank_access(
      const dealii::Utilities::MPI::Partitioner &partitioner)
    {
      // count the number of times a locally owned DoF is accessed by the
      // remote ghost data
      std::vector<unsigned int> touch_count(partitioner.locally_owned_size());
      for (const auto &p : partitioner.import_indices())
        for (unsigned int i = p.first; i < p.second; ++i)
          touch_count[i]++;

      // category 0: DoFs never touched by ghosts
      std::vector<std::vector<unsigned int>> result(1);
      for (unsigned int i = 0; i < touch_count.size(); ++i)
        if (touch_count[i] == 0)
          result.back().push_back(i);

      // DoFs with 1 appearance can be simply categorized by their (single)
      // MPI rank, whereas we need to go an extra round for the remaining DoFs
      // by collecting the owning processes by hand
      std::map<unsigned int, std::vector<unsigned int>>
        multiple_ranks_access_dof;
      const std::vector<std::pair<unsigned int, unsigned int>> &import_targets =
        partitioner.import_targets();
      auto it = partitioner.import_indices().begin();
      for (const std::pair<unsigned int, unsigned int> &proc : import_targets)
        {
          result.emplace_back();
          unsigned int count_dofs = 0;
          while (count_dofs < proc.second)
            {
              for (unsigned int i = it->first; i < it->second;
                   ++i, ++count_dofs)
                {
                  if (touch_count[i] == 1)
                    result.back().push_back(i);
                  else
                    multiple_ranks_access_dof[i].push_back(proc.first);
                }
              ++it;
            }
        }
      Assert(it == partitioner.import_indices().end(), ExcInternalError());

      // Now go from the computed map on DoFs to a map on the processes for
      // DoFs with multiple owners, and append the DoFs found this way to our
      // global list
      std::map<std::vector<unsigned int>,
               std::vector<unsigned int>,
               std::function<bool(const std::vector<unsigned int> &,
                                  const std::vector<unsigned int> &)>>
        dofs_by_rank{[](const std::vector<unsigned int> &a,
                        const std::vector<unsigned int> &b) {
          if (a.size() < b.size())
            return true;
          if (a.size() == b.size())
            {
              for (unsigned int i = 0; i < a.size(); ++i)
                if (a[i] < b[i])
                  return true;
                else if (a[i] > b[i])
                  return false;
            }
          return false;
        }};
      for (const auto &entry : multiple_ranks_access_dof)
        dofs_by_rank[entry.second].push_back(entry.first);

      for (const auto &procs : dofs_by_rank)
        result.push_back(procs.second);

      return result;
    }



    // Compute two vectors, the first indicating the best numbers for a
    // MatrixFree::cell_loop and the second the count of how often a DoF is
    // touched by different cell groups, in order to later figure out DoFs
    // with far reach and those with only local influence.
    template <int dim, typename Number, typename VectorizedArrayType>
    std::pair<std::vector<unsigned int>, std::vector<unsigned char>>
    compute_mf_numbering(
      const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
      const unsigned int                                  component)
    {
      const IndexSet &owned_dofs = matrix_free.get_dof_info(component)
                                     .vector_partitioner->locally_owned_range();
      const unsigned int n_comp =
        matrix_free.get_dof_handler(component).get_fe().n_components();
      Assert(
        matrix_free.get_dof_handler(component).get_fe().n_base_elements() == 1,
        ExcNotImplemented());
      const bool is_fe_q = dynamic_cast<const FE_Q_Base<dim> *>(
        &matrix_free.get_dof_handler(component).get_fe().base_element(0));

      const unsigned int fe_degree =
        matrix_free.get_dof_handler(component).get_fe().degree;
      const unsigned int nn = fe_degree - 1;

      // Data structure used further down for the identification of various
      // entities in the hierarchical numbering of unknowns. The first number
      // indicates the offset from which a given object starts its range of
      // numbers in the hierarchical DoF numbering of FE_Q, and the second the
      // number of unknowns per component on that particular component. The
      // numbers are group by the 3^dim possible objects, listed in
      // lexicographic order.
      std::array<std::pair<unsigned int, unsigned int>,
                 dealii::Utilities::pow(3, dim)>
        dofs_on_objects;
      if (dim == 1)
        {
          dofs_on_objects[0] = std::make_pair(0U, 1U);
          dofs_on_objects[1] = std::make_pair(2 * n_comp, nn);
          dofs_on_objects[2] = std::make_pair(n_comp, 1U);
        }
      else if (dim == 2)
        {
          dofs_on_objects[0] = std::make_pair(0U, 1U);
          dofs_on_objects[1] = std::make_pair(n_comp * (4 + 2 * nn), nn);
          dofs_on_objects[2] = std::make_pair(n_comp, 1U);
          dofs_on_objects[3] = std::make_pair(n_comp * 4, nn);
          dofs_on_objects[4] = std::make_pair(n_comp * (4 + 4 * nn), nn * nn);
          dofs_on_objects[5] = std::make_pair(n_comp * (4 + 1 * nn), nn);
          dofs_on_objects[6] = std::make_pair(2 * n_comp, 1U);
          dofs_on_objects[7] = std::make_pair(n_comp * (4 + 3 * nn), nn);
          dofs_on_objects[8] = std::make_pair(3 * n_comp, 1U);
        }
      else if (dim == 3)
        {
          dofs_on_objects[0] = std::make_pair(0U, 1U);
          dofs_on_objects[1] = std::make_pair(n_comp * (8 + 2 * nn), nn);
          dofs_on_objects[2] = std::make_pair(n_comp, 1U);
          dofs_on_objects[3] = std::make_pair(n_comp * 8, nn);
          dofs_on_objects[4] =
            std::make_pair(n_comp * (8 + 12 * nn + 4 * nn * nn), nn * nn);
          dofs_on_objects[5] = std::make_pair(n_comp * (8 + 1 * nn), nn);
          dofs_on_objects[6] = std::make_pair(n_comp * 2, 1U);
          dofs_on_objects[7] = std::make_pair(n_comp * (8 + 3 * nn), nn);
          dofs_on_objects[8] = std::make_pair(n_comp * 3, 1U);
          dofs_on_objects[9] = std::make_pair(n_comp * (8 + 8 * nn), nn);
          dofs_on_objects[10] =
            std::make_pair(n_comp * (8 + 12 * nn + 2 * nn * nn), nn * nn);
          dofs_on_objects[11] = std::make_pair(n_comp * (8 + 9 * nn), nn);
          dofs_on_objects[12] = std::make_pair(n_comp * (8 + 12 * nn), nn * nn);
          dofs_on_objects[13] =
            std::make_pair(n_comp * (8 + 12 * nn + 6 * nn * nn), nn * nn * nn);
          dofs_on_objects[14] =
            std::make_pair(n_comp * (8 + 12 * nn + 1 * nn * nn), nn * nn);
          dofs_on_objects[15] = std::make_pair(n_comp * (8 + 10 * nn), nn);
          dofs_on_objects[16] =
            std::make_pair(n_comp * (8 + 12 * nn + 3 * nn * nn), nn * nn);
          dofs_on_objects[17] = std::make_pair(n_comp * (8 + 11 * nn), nn);
          dofs_on_objects[18] = std::make_pair(n_comp * 4, 1U);
          dofs_on_objects[19] = std::make_pair(n_comp * (8 + 6 * nn), nn);
          dofs_on_objects[20] = std::make_pair(n_comp * 5, 1U);
          dofs_on_objects[21] = std::make_pair(n_comp * (8 + 4 * nn), nn);
          dofs_on_objects[22] =
            std::make_pair(n_comp * (8 + 12 * nn + 5 * nn * nn), nn * nn);
          dofs_on_objects[23] = std::make_pair(n_comp * (8 + 5 * nn), nn);
          dofs_on_objects[24] = std::make_pair(n_comp * 6, 1U);
          dofs_on_objects[25] = std::make_pair(n_comp * (8 + 7 * nn), nn);
          dofs_on_objects[26] = std::make_pair(n_comp * 7, 1U);
        }

      const auto renumber_func = [](const types::global_dof_index dof_index,
                                    const IndexSet               &owned_dofs,
                                    std::vector<unsigned int>    &result,
                                    unsigned int &counter_dof_numbers) {
        const types::global_dof_index local_dof_index =
          owned_dofs.index_within_set(dof_index);
        if (local_dof_index != numbers::invalid_dof_index)
          {
            AssertIndexRange(local_dof_index, result.size());
            if (result[local_dof_index] == numbers::invalid_unsigned_int)
              result[local_dof_index] = counter_dof_numbers++;
          }
      };

      unsigned int                         counter_dof_numbers = 0;
      std::vector<unsigned int>            dofs_extracted;
      std::vector<types::global_dof_index> dof_indices(
        matrix_free.get_dof_handler(component).get_fe().dofs_per_cell);

      // We now define a lambda function that does two things: (a) determine
      // DoF numbers in a way that fit with the order a MatrixFree loop
      // travels through the cells (variable 'dof_numbers_mf_order'), and (b)
      // determine which unknowns are only touched from within a single range
      // of cells and which ones span multiple ranges (variable
      // 'touch_count'). Note that this process is done by calling into
      // MatrixFree::cell_loop, which gives the right level of granularity
      // (when executed in serial) for the scheduled vector operations. Note
      // that we pick the unconstrained indices in the hierarchical order for
      // part (a) as this makes it easy to identify the DoFs on the individual
      // entities, whereas we pick the numbers with constraints eliminated for
      // part (b). For the latter, we keep track of each DoF's interaction
      // with different ranges of cell batches, i.e., call-backs into the
      // operation_on_cell_range() function.
      const unsigned int        n_owned_dofs = owned_dofs.n_elements();
      std::vector<unsigned int> dof_numbers_mf_order(
        n_owned_dofs, dealii::numbers::invalid_unsigned_int);
      std::vector<unsigned int> last_touch_by_cell_batch_range(
        n_owned_dofs, dealii::numbers::invalid_unsigned_int);
      std::vector<unsigned char> touch_count(n_owned_dofs);

      const auto operation_on_cell_range =
        [&](const MatrixFree<dim, Number, VectorizedArrayType> &data,
            unsigned int &,
            const unsigned int &,
            const std::pair<unsigned int, unsigned int> &cell_range) {
          for (unsigned int cell = cell_range.first; cell < cell_range.second;
               ++cell)
            {
              // part (a): assign beneficial numbers
              for (unsigned int v = 0;
                   v < data.n_active_entries_per_cell_batch(cell);
                   ++v)
                {
                  // get the indices for the dofs in cell_batch
                  if (data.get_mg_level() == numbers::invalid_unsigned_int)
                    data.get_cell_iterator(cell, v, component)
                      ->get_dof_indices(dof_indices);
                  else
                    data.get_cell_iterator(cell, v, component)
                      ->get_mg_dof_indices(dof_indices);

                  if (is_fe_q)
                    for (unsigned int a = 0; a < dofs_on_objects.size(); ++a)
                      {
                        const auto &r = dofs_on_objects[a];
                        if (a == 10 || a == 16)
                          // switch order x-z for y faces in 3d to lexicographic
                          // layout
                          for (unsigned int i1 = 0; i1 < nn; ++i1)
                            for (unsigned int i0 = 0; i0 < nn; ++i0)
                              for (unsigned int c = 0; c < n_comp; ++c)
                                renumber_func(
                                  dof_indices[r.first + r.second * c + i1 +
                                              i0 * nn],
                                  owned_dofs,
                                  dof_numbers_mf_order,
                                  counter_dof_numbers);
                        else
                          for (unsigned int i = 0; i < r.second; ++i)
                            for (unsigned int c = 0; c < n_comp; ++c)
                              renumber_func(
                                dof_indices[r.first + r.second * c + i],
                                owned_dofs,
                                dof_numbers_mf_order,
                                counter_dof_numbers);
                      }
                  else
                    for (const types::global_dof_index i : dof_indices)
                      renumber_func(i,
                                    owned_dofs,
                                    dof_numbers_mf_order,
                                    counter_dof_numbers);
                }

              // part (b): increment the touch count of a dof appearing in the
              // current cell batch if it was last touched by another than the
              // present cell batch range (we track them via storing the last
              // cell batch range that touched a particular dof)
              data.get_dof_info(component).get_dof_indices_on_cell_batch(
                dofs_extracted, cell);
              for (const unsigned int dof_index : dofs_extracted)
                if (dof_index < n_owned_dofs &&
                    last_touch_by_cell_batch_range[dof_index] !=
                      cell_range.first)
                  {
                    ++touch_count[dof_index];
                    last_touch_by_cell_batch_range[dof_index] =
                      cell_range.first;
                  }
            }
        };

      // Finally run the matrix-free loop.
      Assert(matrix_free.get_task_info().scheme ==
               dealii::internal::MatrixFreeFunctions::TaskInfo::none,
             ExcNotImplemented("Renumbering only available for non-threaded "
                               "version of MatrixFree::cell_loop"));

      matrix_free.template cell_loop<unsigned int, unsigned int>(
        operation_on_cell_range, counter_dof_numbers, counter_dof_numbers);

      AssertDimension(counter_dof_numbers, n_owned_dofs);

      return std::make_pair(dof_numbers_mf_order, touch_count);
    }

  } // namespace



  template <int dim,
            int spacedim,
            typename Number,
            typename VectorizedArrayType>
  std::vector<types::global_dof_index>
  compute_matrix_free_data_locality(
    const DoFHandler<dim, spacedim>                    &dof_handler,
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free)
  {
    Assert(matrix_free.indices_initialized(),
           ExcMessage("You need to set up indices in MatrixFree "
                      "to be able to compute a renumbering!"));

    // try to locate the `DoFHandler` within the given MatrixFree object.
    unsigned int component = 0;
    for (; component < matrix_free.n_components(); ++component)
      if (&matrix_free.get_dof_handler(component) == &dof_handler)
        break;

    Assert(component < matrix_free.n_components(),
           ExcMessage("Could not locate the given DoFHandler in MatrixFree"));

    // Summary of the algorithm below:
    // (a) compute renumbering of each DoF in the order the corresponding
    //     object appears in the mf loop -> local_numbering
    // (b) determine by how many cell groups (call-back places in the loop) a
    //     dof is touched -> touch_count
    // (c) determine by how many MPI processes a dof is touched
    //     -> dofs_by_rank_access
    // (d) combine both category types of (b) and (c) and list the indices
    //     according to four categories

    const std::vector<std::vector<unsigned int>> dofs_by_rank_access =
      group_dofs_by_rank_access(
        *matrix_free.get_dof_info(component).vector_partitioner);

    const auto &[local_numbering, touch_count] =
      compute_mf_numbering(matrix_free, component);

    const IndexSet &owned_dofs = matrix_free.get_dof_info(component)
                                   .vector_partitioner->locally_owned_range();
    const types::global_dof_index locally_owned_size = owned_dofs.n_elements();
    AssertDimension(locally_owned_size, local_numbering.size());
    AssertDimension(locally_owned_size, touch_count.size());

    // Create a second permutation to group the unknowns into the following
    // four categories, to be eventually composed with 'local_renumbering'
    enum class Category : unsigned char
    {
      single,     // DoFs only referring to a single batch and MPI process
      multiple,   // DoFs referring to multiple batches (long-range)
      mpi_dof,    // DoFs referring to DoFs in exchange with other MPI processes
      constrained // DoFs without any reference (constrained DoFs)
    };
    std::vector<Category> categories(locally_owned_size);
    for (unsigned int i = 0; i < locally_owned_size; ++i)
      switch (touch_count[i])
        {
          case 0:
            categories[local_numbering[i]] = Category::constrained;
            break;
          case 1:
            categories[local_numbering[i]] = Category::single;
            break;
          default:
            categories[local_numbering[i]] = Category::multiple;
            break;
        }
    for (unsigned int chunk = 1; chunk < dofs_by_rank_access.size(); ++chunk)
      for (const auto i : dofs_by_rank_access[chunk])
        categories[local_numbering[i]] = Category::mpi_dof;

    // Assign numbers to the categories in the order listed in the enum, which
    // is done by starting 'counter' for each category at an offset
    std::array<unsigned int, 4> n_entries_per_category{};
    for (const Category i : categories)
      {
        AssertIndexRange(static_cast<int>(i), 4);
        ++n_entries_per_category[static_cast<int>(i)];
      }
    std::array<unsigned int, 4> counters{};
    for (unsigned int i = 1; i < 4; ++i)
      counters[i] = counters[i - 1] + n_entries_per_category[i - 1];
    std::vector<unsigned int> numbering_categories;
    numbering_categories.reserve(locally_owned_size);
    for (const Category category : categories)
      {
        numbering_categories.push_back(counters[static_cast<int>(category)]);
        ++counters[static_cast<int>(category)];
      }

    // The final numbers are given by the composition of the two permutations
    std::vector<dealii::types::global_dof_index> new_global_numbers(
      locally_owned_size);
    for (unsigned int i = 0; i < locally_owned_size; ++i)
      new_global_numbers[i] =
        owned_dofs.nth_index_in_set(numbering_categories[local_numbering[i]]);

    return new_global_numbers;
  }

} // namespace DoFRenumbering



/*-------------- Explicit Instantiations -------------------------------*/
#include "dofs/dof_renumbering.inst"


DEAL_II_NAMESPACE_CLOSE
